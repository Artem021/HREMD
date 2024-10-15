import re, os

# import engines
import structure
import base_utils

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
