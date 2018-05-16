import math
import numpy as np

class Frame():
    pass

class Lattice():
    def __init__(self, zero = np.array([0,0,0], dtype=float)):
        self.grain = 0.1
        self.zero = zero

    def cast(self, coordinates): 
        moved_coordinates = coordinates - self.zero #moving atom coordinates
        bin = np.array([math.floor(x) for x in moved_coordinates/self.grain]) #putting atom into bin
        return bin

class Simulation():
    pass

class Potential():
    pass

class Replica():
    pass

class Chain():
    def __init__(self, cas_list, seq):
        self.seq = seq #protein sequence
        self.cas = cas_list #list of coordinates (bins) of ca

class PDBreader():

    def read(self, pdb_file):
        seq = [] #protein sequence
        ca_positions = []
        min_ca_pos = np.array([math.inf,math.inf,math.inf], dtype=float)


        with open(pdb_file) as f:
            for line in f:
                if line[0:6].strip() == 'ATOM':
                    if line[12:16].strip() == 'CA':
                        seq.append(line[17:20].strip())
                        x = float(line[30:38].strip())
                        min_ca_pos[0] = min(min_ca_pos[0], x)
                        y = float(line[38:46].strip())
                        min_ca_pos[1] = min(min_ca_pos[1], y)
                        z = float(line[46:54].strip())
                        min_ca_pos[2] = min(min_ca_pos[2], z)
                        ca_positions.append(np.array([x,y,z]))
                elif line[0:6].strip() == 'ENDMDL':
                    break

        lattice = Lattice(min_ca_pos)

        ca_bins = [lattice.cast(x) for x in ca_positions]

        chain = Chain(ca_bins, seq)

        return chain, lattice



