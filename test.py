from base import *


def test1(pdb_file):
    simulation = Simulation(0, f=pdb_file)
    data, cas = simulation.read_pdb_file()
    caster = Cast()
    caster.get_min_ca_coords(cas)
    caster.get_shift(cas)


test1('10gs_protein1.pdb')
