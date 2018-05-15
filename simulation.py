class Lattice:

    def __init__(self, pos, g = 0.1)
        
        
        self.grein = g
        self.zeros = [min(pos[:,0]),min(pos[:,1]),min(pos[:,2])]
    
    def cast(self, coords):
    
        return (int, int, int)

    
    
class PDBReader:


    def read(self, file):
        
        seq = []
        pos = []
        o = open(file, 'r')
        for line in o:
            line = line.split()
            if line[0]=='ATOM':
                if line[2]=='CA':
                    seq.append(line[3])
                    pos.append((float(line[6]), float(line[7]), float(line[8])))
        
        lattice = Lattice(pos)
        chain = Chain(
        #po przeczytaniu 
        return chain, frame, lattice
    
    def to_pdb():
