import json, os, subprocess, sys
from structure import REMD, optimizeFrames
import base_utils

with open(f'config.json','r') as parm:
    parms = json.load(parm)

globalParm, REMDSettings, options  = map(parms.get, ('GLOBAL','REMD_Options','LAMMPS_Options'))

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
    sys.exit()

elem = base_utils.getElementsLmp(DataFile)


vars = {
    'N' : MDsteps,
    'S' : Seed,
    'T' : Temperature
}

parm = {
    'dump_modify' : 'DUMPFILE element '+' '.join(elem),
    'improper_style' : 'umbrella',
    'dihedral_style' : 'harmonic'
    } # TODO: add feature to getXyzfromData()


remdSim = REMD(AlphaRange, DataFile, wd, Nmax=Nmax, seed=Seed, T=Temperature, vars=vars, NPTs = NPTs, delayNPT=delayNPT, Ncores = ncores, selectLastStruc=selectLastStruc,options=options, parm=parm)

remdSim.unrestrictedExchange = unrestrictedExchange
remdSim.selectLastStruc = selectLastStruc
remdSim.changeOrder = changeOrder
remdSim.addWorlds = addWorlds
remdSim.delWorlds = delWorlds
remdSim.Pmin = Pmin
remdSim.Pmax = Pmax
remdSim.Ncheck = Ncheck

remdSim.runREMD(Niter)

optimizeFrames(os.path.join(wd,'structures.xyz'), DataFile, maxp=10, cores=ncores, parm = {
# optimizeFrames('/home/artem/LAMMPS_TEST/macro/22.06/structures.xyz', DataFile, parm = {
    'minimize' : '1.0e-4 1.0e-6 5000 1000', 
    'write_dump' : ' all xyz $t modify element '+' '.join(elem),
    'improper_style' : 'umbrella',
    'dihedral_style' : 'harmonic'
    }, options = {'optimize' : True})
print(0)

