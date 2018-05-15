class Simulation:
    def __init__(self, nrep, file):
        pass

    def run(self):
        pass


class Lattice:
    def __init__(self, grain):
        self.grain = grain
        self.zeros = (0, 0, 0)

    def cast(self, coords):
        pass


class Chain:
    def __init__(self, seq, ):
        pass


class PDBReader:
    def __init__(self, lc):
        self.lc = Lattice()

    def read(self, file):
        seq = []
        pos = []
        o = open(file, "r")
        for line in o.readlines():
            line = line.split()
            if line[0]=='ATOM':
                if line[1]=='CA':
                    seq.append(line[3])
                    pos.append((float(line[6]), float(line[7]), float(line[8])))


class Frame:
    pass


class Potential:
    pass


class Replica:
    pass
