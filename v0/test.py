import os, re
alphaSet = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
# alphaSet = [0.0, 1.0]
reverse = True


order = range(len(alphaSet)-1)
if reverse:
    order = range(len(alphaSet)-2,-1,-1)

for i in order:
    w1 = alphaSet[i]
    w2 = alphaSet[i+1]
    print(f'{w1} --> {w2}')


exit()
path = 'src/lammps.data'

def getXyzfromData(path):
    key = 'Atoms'
    nextKey = 'Bonds'
    natoms = 0
    xyz = []
    with open(path,'r') as fi:
        lines = iter(fi.readlines())
        read = False
        while True:
            line = next(lines)
            if nextKey in line or line == '' or read==True:
                raise RuntimeError('Incorrect data file')
            if key in line:
                while True:
                    # row = next(lines).split()
                    row = next(lines)
                    row = row.split()
                    if row != []:
                        try:
                            *_, x, y, z = row[-3:]
                            xyz.append('    '.join([x,y,z]))
                            natoms+=1
                        except:
                            xyz.insert(0,f'structure from {path}')
                            xyz.insert(0,f'{natoms}')
                            return '\n'.join(xyz)
                        

xyz = getXyzfromData(path)

print(xyz)






