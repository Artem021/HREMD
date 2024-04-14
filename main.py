import json, os, subprocess
from structure import REMD

with open(f'config.json','r') as parm:
    parms = json.load(parm)

globalParm, REMDSettings, options  = map(parms.get, ('GLOBAL','REMD_Options','LAMMPS_Options'))

JobName, CalcDir, Seed, Temperature, EnergyUnits, LAMMPS_PATH, DataFile, \
    XyzFile, MDsteps, Niter, AlphaRange = map(globalParm.get, \
        ('JobName', 'CalcDir', 'Seed', 'Temperature', 'EnergyUnits', 'LAMMPS_PATH',\
        'DataFile', 'XyzFile', 'MDsteps', 'Iterations', 'AlphaRange'))

unrestrictedExchange, selectLastStruc, changeOrder, addWorlds, delWorlds, Pmin, \
    Pmax, Ncheck, Nmax = map(REMDSettings.get, ('unrestrictedExchange',\
        'selectLastStruc', 'changeOrder', 'addWorlds', 'delWorlds', 'Pmin', 'Pmax', \
        'Ncheck', 'Nmax'))

wd = os.path.join(CalcDir,JobName)

os.environ['PATH'] = LAMMPS_PATH
if subprocess.call('lmp -help',shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE):
    raise RuntimeError(f'LAMMPS not found in {LAMMPS_PATH}')

vars = {
    'N' : MDsteps,
    'S' : Seed,
    'T' : Temperature
}

remdSim = REMD(AlphaRange, DataFile, wd, Nmax=Nmax, seed=Seed, T=Temperature, vars=vars, options=options)

remdSim.unrestrictedExchange = unrestrictedExchange
remdSim.selectLastStruc = selectLastStruc
remdSim.changeOrder = changeOrder
remdSim.addWorlds = addWorlds
remdSim.delWorlds = delWorlds
remdSim.Pmin = Pmin
remdSim.Pmax = Pmax
remdSim.Ncheck = Ncheck


remdSim.runREMD(Niter)

