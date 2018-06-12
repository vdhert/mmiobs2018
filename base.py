import numpy as np


class Lattice:
    map = None
    initial_data = None
    min_ca_coords = None

    def __init__(self, data, size):
        self.map = np.zeros(size)
        self.initial_data = data
        self.get_min_ca_coords_and_shift()

    def cast(self, data):
        coords = np.array([res for _, res in self.initial_data])
        return np.floor(coords / self.map.size[0])

    def get_min_ca_coords_and_shift(self):
        for ca in self.initial_data:
            if self.min_ca_coords is None:
                self.min_ca_coords = ca[1]
            else:
                for coord in range(3):
                    if self.min_ca_coords[coord] > ca[1][coord]:
                        self.min_ca_coords[coord] = ca[1][coord]

        for _, coords in self.initial_data:
            coords += self.min_ca_coords
        self.initial_data = self.initial_data


class Simulation:

    def __init__(self, n_rep, f):
        self.n_rep = n_rep
        self.n_steps = 0
        self.file = f

    def read_pdb_file(self):
        with open(self.file, 'r+') as f:
            protein = f.readlines()
            i = 0
            cas = []
            while not protein[i].startswith('ATOM'):
                i += 1
            while protein[i].startswith('ATOM'):
                line = protein[i].split()
                amino_acid = (line[3], np.array(line[6:9]).astype(np.float).round(3))
                if line[2] == 'CA':
                    cas.append(amino_acid)
                i += 1
            return cas


class Frame:
    recs = []
    gcs = []
# gcs = recs


class Chain:
    cas = []
    seq = []
    recs = []
