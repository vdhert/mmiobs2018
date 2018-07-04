from Bio.PDB import PDBParser
import numpy as np
import gzip
import os
from resn_dict import Res
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


# produktem pętli jest lista list w postaci ['1-literowa nazwa reszty', pozycja C-1, pozycja C0, pozycja C1,
# pozycja środka geometrycznego łańcucha bocznego]
def res_from_protein(pdb):
    res_list = []
    d = Res()
    d = d.res_dict
    for model in pdb:
        for chain in model:
            c_before = 0
            c = 0
            for res in chain:
                if res.id[0] == ' ':
                    c_before = c
                    try:
                        c = res['CA'].get_coord()
                    except KeyError:
                        print(res)
                    if type(c_before) != int:
                        sidechain = []
                        for atom in res:
                            if atom.get_name() not in ["N", "CA", "O", "C"]:
                                sidechain.append(atom.get_coord())
                        try:
                            cg = sum(sidechain) / len(sidechain)  # liczenie środka geometrycznego
                        except ZeroDivisionError:
                            cg = 0
                        try:
                            res_parameters = [d[res.resname], c_before, c, c, cg, str(res.id)]
                        except KeyError:
                            res_parameters = ['x', c_before, c, c, cg, res.id]
                        if len(res_list) > 0:
                            res_list[-1][3] = c
                        res_list.append(res_parameters)
    return res_list[:-1]  # ucinamy ostatni - brak położenia C+1


# poniższa funkcja oblicza wektor położenia środka geometrycznego po obrocie aa. Zwraca 1-literową nazwę reszty i
# położenie środka dla C0 = 0, 0, 0.
def turn_cg(res):
    # pattern = [np.array([0., 0., 0.]), np.array([2.3, 2.3, 0.]),np.array([5.6, 0.0, 0.0])]
    # od razu przesuwamy bb_coords tak, by C-1 leżał w początku układu współrzędnych

    # druga wersja - przesuwamy tak, by C0 leżało w początku układu, pattern też jest przesunięty
    pattern = [np.array([-2.3, -2.3, 0.]), np.array([0., 0., 0.]), np.array([3.3, -2.3, 0.0])]
    bb_coords = [res[1] - res[2], res[2] - res[2], res[3] - res[1]]

    # U - macierz obrotu
    V, S, Wt = np.linalg.svd(np.dot(np.transpose(bb_coords), pattern))
    reflect = float(str(float(np.linalg.det(V) * np.linalg.det(Wt))))
    if reflect == -1.0:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]
    U = np.dot(V, Wt)

    # teraz przesuwamy i obracamy cg

    cg = res[4] - res[2]
    new_cg = np.dot(cg, U)

    return res[0], new_cg


def create_cg_dict(source='f', path="./bialka"):
    file_list = []
    if source == 'f':
        for f in os.listdir(path):
            if (f.endswith("pdb") or f.endswith("ent")) and not f.startswith('~'):
                file_list.append(path + "/" + f)
    elif source == 'pdb':
        file_list.append(path)

    cg_dict = {}
    for f in file_list:
        parser = PDBParser()
        pdb = parser.get_structure("test", f)
        #print(f)
        res_list = res_from_protein(pdb)
        for res in res_list:
            aa, cg = turn_cg(res)
            if aa not in cg_dict.keys():
                cg_dict[aa] = [cg]
            else:
                cg_dict[aa].append(cg)

    avg_cg_dict = {k: (sum(v)/len(v)) for k, v in cg_dict.items()}
    #return avg_cg_dict
    return cg_dict


if __name__ == '__main__':

    d = create_cg_dict('pdb', "pdb_test.txt")
    #k, cg, U = turn_cg(['A', np.array([-2.3, -2.3, 0.]), np.array([0., 0., 0.]), np.array([-2.3, -2.3, 0.]), np.array([0., 1., 0.])])
    print(d)

#co chcę wyprodukować
