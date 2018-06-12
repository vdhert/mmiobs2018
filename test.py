from base import *
from kabsch import *


def test1(pdb_file):
    simulation = Simulation(0, f=pdb_file)
    cas = simulation.read_pdb_file()
    lattice = Lattice(cas)
    lattice.cast()
    print(lattice.)
    pattern = [
        [0., 0., 0.],
        [2.3, 2.3, 0.],
        [5.6, 0.0, 0.0]]


test1('1fn3.pdb')
