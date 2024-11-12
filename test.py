import re, os
import pandas as pd
import engines
import structure
import base_utils

TEST_OPTIMIZATION = False
datafile = '/home/md/md/HREMD_v1/toy.lmps'
xyz = '/home/md/md/HREMD_v1/toy_relax.xyz'
# E = f(alpha)
dh = 0.01

alphas = [round(i*dh,2) for i in range(int(1/dh)+1)]
ener = []

elem = base_utils.getElementsLmp(datafile)
# struc = structure.Structure(datafile,'','','','','','', xyz=base_utils.readXYZ(xyz,0))
struc = structure.Structure(datafile,'','','','','','', xyz=None)

for alp in alphas:
    sim = engines.Simulation(alp,
                             '/home/md/md/HREMD_v1/single_point_test',
                             datafile, 
                             parm = {'improper_style' : 'umbrella', 'dihedral_style' : 'harmonic','dump_modify' : 'DUMPFILE element '+' '.join(elem)})
    world = (alp,struc,sim)
    ener.append(structure.REMD.calc_single_point(world,struc))

df = pd.DataFrame.from_dict({'alpha' : alphas, 'energy' : ener})
df.to_excel("/home/md/md/HREMD_v1/toy_unrelaxed.xlsx")

print(0)
exit()

log = '/home/md/md/MDtest/publ/opt_md/remd_x20-rework/prob.txt'
with open(log,'r') as fi:
    iterlines = iter(fi.readlines())
patt = 'Iteration'
res = {}
while True:
    try:
        line = next(iterlines)
    except StopIteration:
        break
    if re.search(patt,line):
        group = '_'.join(re.split('; |= ',line)[1::2])
        next(iterlines)
        delta = float(next(iterlines).split()[-1])
        if res.get(group):
            res[group].append(delta)
        else:
            res[group] = [delta]

for i in res:
    print(len(res[i]))
df = pd.DataFrame.from_dict(res)
df.to_excel("prob.xlsx")
print(0)
exit()

if TEST_OPTIMIZATION:
    # sim = engines.Simulation(1.0,'/home/md/md/HREMD_v1/test/','/home/md/md/HREMD_v1/test/uniform_polymer_soft.lmps','/home/md/md/HREMD_v1/test/structures-692-check.xyz',options={'optimize':True})
    # sim._writeInput()
    ncores = 24
    wd = '/home/md/md/PES-185_dens/16.07-250_2500-Mn-T600/'
    DataFile = '/home/md/md/PES-185_dens/16.07-250_2500-Mn-T600/uniform_polymer_soft.lmps'
    elem = base_utils.getElementsLmp(DataFile)



    structure.optimizeFrames(os.path.join(wd,'structures.xyz'), DataFile, maxp=10, cores=ncores, parm = {
    # optimizeFrames('/home/artem/LAMMPS_TEST/macro/22.06/structures.xyz', DataFile, parm = {
        'minimize' : '1.0e-4 1.0e-6 5000 1000', 
        'write_dump' : ' all xyz $t modify element '+' '.join(elem),
        'improper_style' : 'umbrella',
        'dihedral_style' : 'harmonic'
        }, options = {'optimize' : True})

    print(0)
