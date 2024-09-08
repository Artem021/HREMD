import os, time, sys, shutil
import subprocess
from multiprocessing import Pool

# root directory
WD = '/home/md/md/PES-185_dens/08.08-mono_8_68/Tg/'
# WD = '/home/md/md/PAM/Tg'
try:
    os.mkdir(WD)
except:
    pass
# datafile with relaxed cell
initialData = '/home/md/md/PES-185_dens/08.08-mono_8_68/Tg/PES-185-initial.data'
# initialData = '/home/md/md/PAM/PAM-initial.data'

# pressure range
# p = [26000, 32000, 38000, 44000]
p = [1, 100, 1000, 10000]
# starting temperature
T1 = 650 # 650
# temperature after cooling
T2 = 150 # 150

# production steps
nprod = 20000 # 20000
nnpt = 50000 # 50000
ncool = 500000 # 500000

# repeat production n times
groups = 5

# mpi tasks per calculation
CORES = 12

# maximum threads
MAXPROC = 20

# lmps exe
LMP = 'lmp'


# npt_tmpl = '''units real
# dimension      3
# processors     * * *
# boundary       p p p
# atom_style     full

# # pair_style      lj/cut/coul/long 10
# pair_style lj/cut/coul/long/soft 1.0 0.5 10.0 8.0 8.0
# bond_style      harmonic
# angle_style     harmonic
# # dihedral_style  harmonic
# # improper_style umbrella
# dihedral_style opls
# improper_style cvff
# # kspace_style    pppm 1.0e-4

# pair_modify     tail yes mix arithmetic

# # read the original positions for the atoms of the molecule, as well as the same parameter file as previously
# read_data %(datafile)s# nvt600.data
# kspace_style    pppm 1.0e-4

# fix      1 all npt temp %(T)s %(T)s 100 iso %(P)s %(P)s 1000
# dump mydmp all atom 1000 npt600_%(P)s.lammpstrj
# thermo 1000
# thermo_style   custom step temp pe press vol lx density
# variable mytemp equal temp
# fix myat1 all ave/time 10 10 200 v_mytemp file temperature.dat

# timestep 1
# run            %(N)s
# write_data npt600_%(P)s.data

# unfix          1'''

# prod_templ = '''units real
# dimension      3
# processors     * * *
# boundary       p p p
# atom_style     full
# # pair_style      lj/cut/coul/long 10
# pair_style lj/cut/coul/long/soft 1.0 0.5 10.0 8.0 8.0
# bond_style      harmonic
# angle_style     harmonic
# # dihedral_style  harmonic
# # improper_style umbrella
# dihedral_style opls
# improper_style cvff
# # kspace_style    pppm 1.0e-4
# neigh_modify one 10000
# pair_modify     tail yes mix arithmetic

# # read the original positions for the atoms of the PE molecule, as well as the same parameter file as previously
# read_data %(datafile)s # npt_1_32000.data
# kspace_style    pppm 1.0e-4

# fix      1 all npt temp %(temp)s %(temp)s 100 iso %(press)s %(press)s 1000 # p = 1; t = 650
# #dump mydmp all atom 1000 prod.lammpstrj
# thermo 1000
# thermo_style   custom step temp pe press vol lx density
# variable mytemp equal temp
# fix myat1 all ave/time 10 10 200 v_mytemp file temperature.dat

# timestep 1
# %(runstring)s
# unfix          1'''

# cool_templ = '''#package gpu 1

# units real
# dimension      3
# processors     * * *
# boundary       p p p
# atom_style     full
# # pair_style      lj/cut/coul/long/gpu 10
# # pair_style      lj/cut/coul/long 10
# pair_style lj/cut/coul/long/soft 1.0 0.5 10.0 8.0 8.0
# bond_style      harmonic
# angle_style     harmonic
# dihedral_style opls
# improper_style cvff
# # dihedral_style  harmonic
# # improper_style umbrella
# # kspace_style    pppm 1.0e-4

# #special_bonds   dreiding
# pair_modify     tail yes mix arithmetic

# # read the original positions for the atoms of the PEG molecule, as well as the same parameter file as previously
# read_data %(datafile)s # 8000_3.data
# kspace_style    pppm 1.0e-4

# # fix      1 all npt temp %(T1)s %(T2)s 100 iso %(P)s %(P)s 1000
# fix      1 all npt temp %(T1)s %(T2)s 100 iso 0 0 1000
# dump mydmp all atom 1000 Cooling.lammpstrj
# thermo 1000
# thermo_style   custom step temp pe press vol lx density
# variable mytemp equal temp
# fix myat1 all ave/time 10 10 200 v_mytemp file temperature.dat

# timestep 1
# run            %(N)s # 500000
# write_data cooling_%(datafile)s # cooling_8000_3.data

# unfix          1'''


npt_tmpl = '''units real
dimension      3
processors     * * *
boundary       p p p
atom_style     full

# pair_style      lj/cut/coul/long 10
pair_style lj/cut/coul/long/soft 1.0 0.5 10.0 8.0 8.0
bond_style      harmonic
angle_style     harmonic
dihedral_style  harmonic
improper_style umbrella
kspace_style    pppm 1.0e-4

pair_modify     tail yes mix arithmetic

# read the original positions for the atoms of the molecule, as well as the same parameter file as previously
read_data %(datafile)s# nvt600.data

fix      1 all npt temp %(T)s %(T)s 100 iso %(P)s %(P)s 1000
dump mydmp all atom 1000 npt600_%(P)s.lammpstrj
thermo 1000
thermo_style   custom step temp pe press vol lx density
variable mytemp equal temp
fix myat1 all ave/time 10 10 200 v_mytemp file temperature.dat

timestep 1
run            %(N)s
write_data npt600_%(P)s.data

unfix          1'''

prod_templ = '''units real
dimension      3
processors     * * *
boundary       p p p
atom_style     full
# pair_style      lj/cut/coul/long 10
pair_style lj/cut/coul/long/soft 1.0 0.5 10.0 8.0 8.0
bond_style      harmonic
angle_style     harmonic
dihedral_style  harmonic
improper_style umbrella
kspace_style    pppm 1.0e-4
neigh_modify one 10000
pair_modify     tail yes mix arithmetic

# read the original positions for the atoms of the PE molecule, as well as the same parameter file as previously
read_data %(datafile)s # npt_1_32000.data

fix      1 all npt temp %(temp)s %(temp)s 100 iso %(press)s %(press)s 1000 # p = 1; t = 650
# dump mydmp all atom 1000 prod.lammpstrj
thermo 1000
thermo_style   custom step temp pe press vol lx density
variable mytemp equal temp
fix myat1 all ave/time 10 10 200 v_mytemp file temperature.dat

timestep 1
%(runstring)s
unfix          1'''

cool_templ = '''#package gpu 1

units real
dimension      3
processors     * * *
boundary       p p p
atom_style     full
# pair_style      lj/cut/coul/long/gpu 10
# pair_style      lj/cut/coul/long 10
pair_style lj/cut/coul/long/soft 1.0 0.5 10.0 8.0 8.0
bond_style      harmonic
angle_style     harmonic
dihedral_style  harmonic
improper_style umbrella
kspace_style    pppm 1.0e-4

# special_bonds   dreiding
pair_modify     tail yes mix arithmetic

# read the original positions for the atoms of the PEG molecule, as well as the same parameter file as previously
read_data %(datafile)s # 8000_3.data

# fix      1 all npt temp %(T1)s %(T2)s 100 iso %(P)s %(P)s 1000
fix      1 all npt temp %(T1)s %(T2)s 100 iso 0 0 1000
dump mydmp all atom 1000 Cooling.lammpstrj
thermo 1000
thermo_style   custom step temp pe press vol lx density
variable mytemp equal temp
fix myat1 all ave/time 10 10 200 v_mytemp file temperature.dat

timestep 1
run            %(N)s # 500000
write_data cooling_%(datafile)s # cooling_8000_3.data

unfix          1'''


def runLmp(lmp, inpfile):
    wd = os.path.dirname(inpfile)
    os.chdir(wd)
    return os.system('mpirun -n %s %s -i %s' % (CORES, lmp, inpfile))

def mkdir(wd):
    try:
        os.mkdir(wd)
    except:
        shutil.rmtree(wd,ignore_errors=True)
        os.mkdir(wd)
    print(wd)

def runNPT(wd,data,template,nstep,rangeP,T,from_scratch=False):
    results = {}
    assert os.path.exists(data), 'datafile not found, could not run NPT'
    exc = 0
    save = False
    if os.path.exists(os.path.join(wd,'success')) and not from_scratch:
        print('Avoid overwriting a previous successful calculation')
        save = True
    else:
        print('Running calculation from scratch')
        mkdir(wd)
    os.chdir(wd)
    args = []
    base = os.path.basename(os.getcwd())
    for p in rangeP:
        subdir = '%s_%s' % (base,p)
        if not save:
            mkdir(subdir)
            os.chdir(subdir)
            cont = template % {'datafile' : os.path.basename(data),  'T':T, 'P': p, 'N':nstep}
            with open('input_npt600.lammps','w') as fo:
                fo.write(cont)
            shutil.copy2(data,'.')
            args.append((LMP, os.path.abspath('input_npt600.lammps')))
            os.chdir('..')
        results[p] = os.path.abspath(os.path.join(subdir,'npt600_%s.data' % p))
    if save:
        return 0,results
    with Pool(MAXPROC) as pool:
        res = pool.starmap(runLmp, args)
    if any(res):
        exc = 1
    if all(res):
        exc = 2
    if exc==0:
        os.chdir(wd)
        open('success','w').close()
    return exc,results

def runProd(wd, data, template, ngroups, nstep, rangeP, T, from_scratch=False):
    results = {}
    exc = 0
    save = False
    if os.path.exists(os.path.join(wd,'success')) and not from_scratch:
        print('Avoid overwriting a previous successful calculation')
        save = True
    else:
        print('Running calculation from scratch')
        mkdir(wd)
    os.chdir(wd)
    args = []
    for p in rangeP:
        subdir = 'Prod_%s' % p
        if not save:
            df = data[p]
            if not os.path.exists(df):
                print('Missing data for p = %s atm, skipping' % p)
                continue
            mkdir(subdir)
            os.chdir(subdir)
            runstring = ''.join(['run            %(nstep)s\nwrite_data %(p)s_%(i)s.data\n\n' % {'nstep':nstep,'p':p, 'i':i } for i in range(1, ngroups+1)])
            cont = template % {'datafile' : os.path.basename(df), 'runstring' : runstring, 'temp':T, 'press': 1}
            with open('input_prod.lammps','w') as fo:
                fo.write(cont)
            shutil.copy2(df,'.')
            args.append((LMP, os.path.abspath('input_prod.lammps')))
            os.chdir('..')
        results[p] = [os.path.abspath(os.path.join(subdir,'%s_%s.data' % (p,i))) for i in range(1,ngroups+1)]
    if save:
        return 0, results
    with Pool(MAXPROC) as pool:
        res = pool.starmap(runLmp, args)
    if any(res):
        exc = 1
    if all(res):
        exc = 2
    if exc==0:
        os.chdir(wd)
        open('success','w').close()
    return exc,results

def runCool(wd, data, template, ngroups, nstep, rangeP, T1, T2, eigenP=False):
    # assert os.path.exists(data), 'datafile not found, production run failed'
    exc = 0
    mkdir(wd)
    os.chdir(wd)
    args = []
    for p in rangeP:
        try:
            datafiles = data[p]
        except:
            print('No data for pressure = %s atm' % p)
            continue
        _p = 0
        if eigenP:
            _p = p
        # datafiles = data.get(p,[])
        if not all([os.path.exists(file) for file in datafiles]) and len(datafiles) != 0:
            print('Missing data for pressure=%s atm, skip' % p)
            continue
        for n in range(1,ngroups+1):
            subdir = 'Cool_%s_%s' % (p,n)
            mkdir(subdir)
            os.chdir(subdir)
            cont = template % {'datafile' : os.path.basename(datafiles[n-1]),  'T1':T1, 'T2':T2, 'P': _p, 'N':nstep}
            with open('input_cooling.lammps','w') as fo:
                fo.write(cont)
            shutil.copy2(datafiles[n-1],'.')
            args.append((LMP, os.path.abspath('input_cooling.lammps')))
            os.chdir('..')
    with Pool(MAXPROC) as pool:
        res = pool.starmap(runLmp, args)
    if any(res):
        exc = 1
    if all(res):
        exc = 2
    return exc
    
# relax cell p=1 atm --> npt under p =[1,100,1000,10000], t1 --> production p=1, t1 --> cooling p, [t1->t2]


# 0. run preopt, make initial data file

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# + main script for Glass Transition Temperature determination +
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# 1. NPT under different pressures
t0 = time.time()
rc,datafiles = runNPT(os.path.join(WD,'./npt600'),initialData,npt_tmpl,nnpt,p,T1)
t1 = time.time()
if rc==0:
    t = (t1-t0)/60
    print('NPT done successfully in %.0f min\n' % t)
elif rc==2:
    print('NPT failed for some reason, exiting\n')
    sys.exit()
else:
    print('some calculations failed, try to continue\n')

# 2. PRODUCTION

t0 = time.time()
rc,datafiles = runProd(os.path.join(WD,'./Production'), datafiles, prod_templ, groups, nprod, p, T1)
t1 = time.time()
if rc==0:
    t = (t1-t0)/60
    print('production done successfully in %.0f min\n' % t)
elif rc==2:
    print('Production failed for some reason, exiting')
    sys.exit()
else:
    print('some calculations failed, try to continue\n')

# 3. Cooling
# print(datafiles)
t0 = time.time()

# dtmp = {1: ['/home/artem/LAMMPS/Tglass/test_pes210/Production/Prod_1/1_1.data', '/home/artem/LAMMPS/Tglass/test_pes210/Production/Prod_1/1_2.data', '/home/artem/LAMMPS/Tglass/test_pes210/Production/Prod_1/1_3.data', '/home/artem/LAMMPS/Tglass/test_pes210/Production/Prod_1/1_4.data', '/home/artem/LAMMPS/Tglass/test_pes210/Production/Prod_1/1_5.data'], 100: ['/home/artem/LAMMPS/Tglass/test_pes210/Production/Prod_100/100_1.data', '/home/artem/LAMMPS/Tglass/test_pes210/Production/Prod_100/100_2.data', '/home/artem/LAMMPS/Tglass/test_pes210/Production/Prod_100/100_3.data', '/home/artem/LAMMPS/Tglass/test_pes210/Production/Prod_100/100_4.data', '/home/artem/LAMMPS/Tglass/test_pes210/Production/Prod_100/100_5.data'], 1000: ['/home/artem/LAMMPS/Tglass/test_pes210/Production/Prod_1000/1000_1.data', '/home/artem/LAMMPS/Tglass/test_pes210/Production/Prod_1000/1000_2.data', '/home/artem/LAMMPS/Tglass/test_pes210/Production/Prod_1000/1000_3.data', '/home/artem/LAMMPS/Tglass/test_pes210/Production/Prod_1000/1000_4.data', '/home/artem/LAMMPS/Tglass/test_pes210/Production/Prod_1000/1000_5.data'], 10000: ['/home/artem/LAMMPS/Tglass/test_pes210/Production/Prod_10000/10000_1.data', '/home/artem/LAMMPS/Tglass/test_pes210/Production/Prod_10000/10000_2.data', '/home/artem/LAMMPS/Tglass/test_pes210/Production/Prod_10000/10000_3.data', '/home/artem/LAMMPS/Tglass/test_pes210/Production/Prod_10000/10000_4.data', '/home/artem/LAMMPS/Tglass/test_pes210/Production/Prod_10000/10000_5.data']}
rc = runCool(os.path.join(WD,'./Cooling'), datafiles, cool_templ, groups, ncool, p, T1, T2, eigenP=True)
t1 = time.time()
if rc==0:
    t = (t1-t0)/60
    print('Cooling done successfully in %.0f min\n' % t)
elif rc==2:
    print('Some cooling jobs failed')
    sys.exit()
else:
    print('Cooling error\n')

print(rc)
