import numpy as np


class Cast:
    map = None
    min_ca_coords = None

    def __init__(self):
        self.map = np.zeros((100, 100, 100))

    def create_cast(self, data):
        return

    def get_min_ca_coords(self, data):
        for ca in data:
            if self.min_ca_coords is None:
                self.min_ca_coords = ca[1]
            else:
                for coord in range(3):
                    if self.min_ca_coords[coord] > ca[1][coord]:
                        self.min_ca_coords[coord] = ca[1][coord]

    def get_shift(self, data):
        pass  # TODO
        # for ca in data:
        #     print(ca[1], ca[1] - self.min_ca)


class Simulation:

    def __init__(self, n_rep, f):
        self.n_rep = n_rep
        self.n_steps = 0
        self.file = f

    # Caster chce dostac liste wszystkich wegli alfa
    # i to on ma wyznaczyc minimum i dzielic na biny

    def read_pdb_file(self):
        with open(self.file, 'r+') as f:
            protein = f.readlines()[2:]
            i = 0
            data, cas = [], []
            while protein[i].startswith('ATOM'):
                line = protein[i].split()
                amino_acid = (line[3], np.array(line[6:9]).astype(np.float).round(1))
                if line[2] == 'CA':
                    cas.append(amino_acid)
                else:
                    data.append(amino_acid)
                i += 1
            return data, cas



