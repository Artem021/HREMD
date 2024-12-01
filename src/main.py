
import json, os, subprocess, sys
from structure import REMD, optimizeFrames
import base_utils
import argparse
# reference lammps 
# /home/users/artem_k/lammps-static/bin/

# apptainer xrd
# /opt/lammps/build/

parser = argparse.ArgumentParser(description='aREMD with XRD metadynamics support')
parser.add_argument('input', type=str,
                    help='Location of json file with input options')
args = parser.parse_args()

print(f'Running REMD simulation from {args.input}')
with open(args.input,'r') as parm:
    parms = json.load(parm)



# with open(f'config.json','r') as parm:
#     parms = json.load(parm)
#TODO: pass as kwargs to REMD()
globalParm, REMDSettings, options, parm_lmp, plumed  = map(parms.get, ('GLOBAL','REMD_Options','LAMMPS_Options','LAMMPS_Parameters','PLUMED_Options'))

JobName, CalcDir, Seed, Temperature, EnergyUnits, LAMMPS_PATH, DataFile, \
    XyzFile, MDsteps, Niter, AlphaRange = map(globalParm.get, \
        ('JobName', 'CalcDir', 'Seed', 'Temperature', 'EnergyUnits', 'LAMMPS_PATH',\
        'DataFile', 'XyzFile', 'MDsteps', 'Iterations', 'AlphaRange'))

unrestrictedExchange, selectLastStruc, changeOrder, addWorlds, delWorlds, Pmin, \
    Pmax, Ncheck, Nmax, NPTs, delayNPT, ncores = map(REMDSettings.get, ('unrestrictedExchange',\
        'selectLastStruc', 'changeOrder', 'addWorlds', 'delWorlds', 'Pmin', 'Pmax', \
        'Ncheck', 'Nmax', 'withNPT', 'delayNPT', 'cores_per_replica'))

wd = os.path.join(CalcDir,JobName)

os.environ["PATH"] += os.pathsep + LAMMPS_PATH
try:
    subprocess.call(['lmp','-help'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
except FileNotFoundError:
    print(f'LAMMPS executable not found, check path: {LAMMPS_PATH}')
    # sys.exit()

elem = base_utils.getElementsLmp(DataFile)

vars = {
    'N' : MDsteps,
    'S' : Seed,
    'T' : Temperature
}

parm = {
    'dump_modify' : 'DUMPFILE element '+' '.join(elem)
#    'improper_style' : 'umbrella',
#    'dihedral_style' : 'harmonic'
    } # TODO: add feature to getXyzfromData()
parm.update(parm_lmp)

remdSim = REMD(AlphaRange, DataFile, wd, Nmax=Nmax, seed=Seed, T=Temperature, vars=vars, NPTs = NPTs, delayNPT=delayNPT, Ncores = ncores, selectLastStruc=selectLastStruc,options=options, parm=parm,plumed=plumed)

remdSim.unrestrictedExchange = unrestrictedExchange
remdSim.selectLastStruc = selectLastStruc
remdSim.changeOrder = changeOrder
remdSim.addWorlds = addWorlds
remdSim.delWorlds = delWorlds
remdSim.Pmin = Pmin
remdSim.Pmax = Pmax
remdSim.Ncheck = Ncheck

remdSim.runREMD(Niter)

_parm = parm.copy()
_parm.update({
    'minimize' : '1.0e-4 1.0e-6 5000 1000', 
    'write_dump' : ' all xyz $t modify element '+' '.join(elem)})
_options = options.copy()
_options.update({'optimize' : True})

MINIMIZE_FRAMES = False
if MINIMIZE_FRAMES:
    optimizeFrames(os.path.join(wd,'structures.xyz'), DataFile, maxp=10, cores=ncores, parm = _parm, options = _options)
    optimizeFrames(os.path.join(wd,'structures_w1.xyz'), DataFile, maxp=10, cores=ncores, parm = _parm, options = _options)
else:
    print('Optimization of frames omitted')
print(0)

