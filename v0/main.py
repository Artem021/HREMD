import os,math,random,shutil,time,queue,json
from multiprocessing import Process,Queue, Pool
import thread
import results as postprocess
import powder


calc_dxrd = False


#---------- read parameters ----------#
assert os.path.isfile('params.json'), 'Check name of parameter file!'
with open(f'params.json','r') as parm:
    parms = json.load(parm)

general_parm, MD_parm, aREMD_parm  = \
    map(parms.get, ('GENERAL','MD','aREMD'))

job_name, topology_file, project_name, xyz_start, engine, seed  = \
    map(general_parm.get, ('JOB_NAME', 'TOP','PROJECT_NAME','XYZ', 'engine', 'seed'))

alpha_set, xyz_found, N_max, exchange_all_pairs, select_last_frame,\
    BETA, EXCHANGE_STAT, ENERGY_STAT = \
    map(aREMD_parm.get, ('alpha_set', 'XYZ_SHARED','NMAX','GREEDY_EXCHANGE',\
                         'USE_LAST_STRUC','BETA','EXCHANGE_STAT','ENERGY_STAT'))
MD_steps, T, cell = \
    map(MD_parm.get, ('md_steps', 'TEMPERATURE','CELL'))


BETA = 1/(273.15 * 1.987204259 * 10**-3) # lammps return E in kcal/mol

#----------set directory----------#
assert len(alpha_set) > 2, 'Invalid number of worlds!'

Nxyz = 0
wd = '../' + project_name
rinp = None

try:
    os.mkdir(wd)
except:
    pass
if engine !='lammps':
    shutil.copyfile(os.path.join('src', xyz_start), os.path.join(wd, xyz_start))
else:
    xyz_start = None
shutil.copyfile(os.path.join('src', topology_file), os.path.join(wd, topology_file))
if engine == 'cp2k':
    inpfile = 'cp2k.inp'
    trajname = 'PAM-pos-1.xyz'
elif engine == 'xtb':
    inpfile = 'xtb.inp'
    trajname = 'xtb.trj'
elif engine == 'lammps':
    inpfile = 'lammps.inp'
    trajname = 'traj.xyz'
    rinp = 'lammps.restart'
    rname = 'binary.restart'
    shutil.copyfile(os.path.join('src', rinp), os.path.join(wd, rinp))
shutil.copyfile(os.path.join('src', inpfile), os.path.join(wd, inpfile))
os.chdir(wd)
try:
    os.mkdir(job_name)
except:
    shutil.rmtree(job_name, ignore_errors=True)
    os.mkdir(job_name)
try:
    shutil.copyfile(xyz_found,os.path.join(job_name,xyz_found))
    with open(xyz_found,'r') as fi: # <-- N of initial structures
        natoms = int(fi.readline())
        nlines = sum(1 for i in fi if i != '\n') + 1
        Nxyz = int(nlines/(2+natoms))
    print(f'{Nxyz} initial structures were read')
except:
    Nxyz = 0
    print('initial structures were not provided')

if engine !='lammps':
    shutil.copyfile(xyz_start,os.path.join(job_name,xyz_start))
shutil.copyfile(topology_file,os.path.join(job_name,topology_file))
shutil.copyfile(inpfile,os.path.join(job_name,inpfile))
if engine =='lammps':
    shutil.copyfile(rinp,os.path.join(job_name,rinp))
os.chdir(job_name)

#----------create Worlds----------#
worlds = []
for alp in alpha_set:
    world = thread.World(
        BASE_DIR = os.path.abspath(os.getcwd()),
        alpha = alp,
        TOP = topology_file,
        XYZ = xyz_start,
        INP = inpfile,
        TRJ = trajname,
        RST = rinp,
        RSTBIN = rname,
        XYZF=xyz_found,
        select_last_frame = select_last_frame,
        engine = engine
    )
    worlds.append(world)

# start of aREMD loop
exchange_set = ['{:.3e}'.format(i) for i in alpha_set]
exchange_ord = {alp:exchange_set.index(alp) for alp in exchange_set}
with open(EXCHANGE_STAT,'w') as exc_stat:
    exc_stat.write('   '.join([str(exchange_ord[a]) for a in exchange_set]) + '\n')
with open(ENERGY_STAT,'w') as e_stat:
    e_stat.write('Type'.center(25) + ''.join([a.center(25) for a in exchange_set]) + '\n')
comment = ''
with open('comments.txt','w') as comm:
    comm.write(comment)

world1 = worlds[-1]
N_i = 0
tcalc = 0
t0 = time.time()
random.seed(seed)
while N_i < N_max:
    # 1. parallel launch of MDs
    procs = []
    main_queue = Queue()
    results = []
    do_restart = N_i > 0 # TODO move inside World
    for i,w in enumerate(worlds):
        # if w.shake:
        #     w.shake_cell(N_i)
        #     # w.set_pot()
        calc = thread.MD_calc(w.DIR,w.XYZ,w.ARGS,w.CMD_LINE,w.INP,w.OUT,w.ERR,w.TRJ,w.RST,do_restart,w.alpha,Nxyz,w.XYZ_FOUND,w.engine)
        proc = Process(target = calc.start_job,args=(main_queue,i))
        # proc = Process(target = thread2.start_XTB,args = (w.XYZ,w.ARGS,w.CMD_LINE,w.INP,w.OUT,w.ERR,w.TRJ,w.alpha,main_queue,i))
        procs.append(proc)
        proc.start()
    # loop to prevent deadlock, explained here:
    # stackoverflow.com/questions/31665328/
    stop_wait = False
    while not stop_wait:
        try:
            result = main_queue.get(False, 0.01)
            results.append(result)
        except queue.Empty:
            pass
        all_done = all([p.exitcode!=None for p in procs])
        stop_wait = all_done & main_queue.empty()
    results = sorted(results, key=lambda res: res['index'])
    for i,w in enumerate(worlds):
        w.set_results(results[i])
    tcalc+=max([r['time']for r in results])
    # if all MDs failed, we exit
    if all([w.error for w in worlds]):
        print(f'All calculations failed, exiting (cycle {N_i})')
        break
    # 2. Metropolis–Hastings part of algorithm
    print(f'Cycle {N_i}')
    for w in worlds:
        w.exchanged = False
    for i in range(len(worlds)-1):
        w1 = worlds[i]
        w2 = worlds[i+1]
        if (w1.error | w2.error):
            continue
        if w1.exchanged: # if we want distant neighbors to have the probability of being exchanged
            if exchange_all_pairs:
                struc1 = w1.get_structure_after()  
            else:
                continue
        else:
            struc1 = w1.get_structure_before()
        if w2.exchanged:
            if exchange_all_pairs:
                struc2 = w2.get_structure_after()  
            else:
                continue
        else:
            struc2 = w2.get_structure_before()
        alp1 = w1.alpha
        alp2 = w2.alpha
        e1 = struc1['Energy']
        e2 = struc2['Energy']
        rel_e = -e1*alp1 -e2*alp2 + e2*alp1 + e1*alp2 # delta = (e1-e2)*(alp2-alp1)
        delta = max(-400,min(400,rel_e * BETA)) # returns -400 if rel_e==nan
        if math.isnan(delta):
            continue
        p = math.exp(-delta)
        p0 = random.random()
        exchange = (p>=p0)
        # print(f'CYCLE {N}') # remove
        # print('=BEFORE EXCHANGE=')
        # print(f'exchange attempt: {alp1} and {alp2}')
        # print(f'Erel = {rel_e} a.u.; delta = {delta}')
        # print(f'p = {p}; p0 = {p0}')
        # print(f'exchange: {exchange}')

        print(f'exchange attempt: {alp1} and {alp2}')
        print(f'Etot({alp1}) = {e1} a.u. and Etot({alp2}) = {e2} a.u.')
        print(f'Index({alp1}) = {struc1["index"]} and index({alp2}) = {struc2["index"]}')
        print(f'Relative energy and Delta: {rel_e} and {delta}')
        print(f'p and p0: {p} and {p0}')
        print(f'exchange: {exchange}')
        if exchange:
            w1.set_structure_after(struc2)
            w2.set_structure_after(struc1)
            if engine != 'lammps':
                w1.update_xyz()
                w2.update_xyz()
            else:
                r1 = w1.get_binrestart_name()
                r2 = w2.get_binrestart_name()
                shutil.copy(r1,r1+f'-{N_i}')
                shutil.copy(r2,r2+f'-{N_i}')
                shutil.copy(r1,r2)
                shutil.copy(r2+f'-{N_i}',r1)
                
            w1.exchanged = True
            w2.exchanged = True
            # collect exchange statistics
            i1 = exchange_set.index('{:.3e}'.format(alp1))
            i2 = exchange_set.index('{:.3e}'.format(alp2))
            atmp = exchange_set[i1]
            exchange_set[i1] = exchange_set[i2]
            exchange_set[i2] = atmp
        else:
            w1.set_structure_after(struc1)
            w2.set_structure_after(struc2)
            if engine != 'lammps':
                w1.update_xyz()
                w2.update_xyz()
            else:
                r1 = w1.get_binrestart_name()
                r2 = w2.get_binrestart_name()
                shutil.copy(r1,r1+f'-{N_i}')
                shutil.copy(r2,r2+f'-{N_i}')
                shutil.copy(r1,r2)
                shutil.copy(r2+f'-{N_i}',r1)

    with open(EXCHANGE_STAT,'a') as exc_stat:
        exc_stat.write('   '.join([str(exchange_ord[a]) for a in exchange_set]) + '\n')
    with open(ENERGY_STAT,'a') as e_stat:
        e_stat.write('ENERGY'.center(25) + ''.join(['{:.8e}'.format(w.struc_before_exchange['Energy']).center(25) for w in worlds]) + '\n')
        e_stat.write('VBIAS'.center(25) + ''.join(['{:.8e}'.format(w.struc_before_exchange['Vbias']).center(25) for w in worlds]) + '\n')
    # 3. update shared XYZ
    world1.update_shared_xyz()
    Nxyz+=1
    
    # 4. prepare for next step
    for w in worlds:
        if w.error:
            #TODO: update .xyz with previos/from another world
            continue
        w.rename_trj_file()
    N_i+=1

    # 5. change order of exchange attempts
    # worlds.reverse() # TODO: better understand this feature

# postprocessing
# if RMSD_REF != None:
#     postprocess.calc_rmsd(RMSD_REF,os.path.abspath(os.getcwd()),'rmsd.txt')
#     postprocess.calc_rmsd(RMSD_REF,os.path.abspath(os.getcwd()),'dihedral.txt',do_rmsd=False)



if calc_dxrd:
    CELL={
                'a':11.8050,
                'b':17.1640,
                'c':7.3930,
                'alpha':90.00,
                'beta':90.00,
                'gamma':90.00,
                'SpGr': 'Pcab'
            }
    MAX_PROC = 12
    print('HREMD is done, now processing trajectory:\n')

    print('1. Calculating XRD difference:')

    txrd=0

    txrd+=powder.calc_xrd('structures.xyz', CELL, 'XRD-DIFF.txt', '/home/artem/HREMD/src/paracetamol-orthorhombic-P1.cif', MAX_PROC)
    txrd+=powder.calc_xrd('structures.xyz', CELL, 'XRD-DIFF-SELF.txt', os.path.abspath('structures-1.cif'), MAX_PROC)

    print('2. Start geometry optimization:')
    topt=powder.optimize_all('structures.xyz', MAX_PROC)

    os.chdir('opt')
    txrd+=powder.calc_xrd('structures-opt.xyz', CELL, 'XRD-DIFF-OPT.txt', '/home/artem/HREMD/src/paracetamol-orthorhombic-P1.cif', MAX_PROC)
    txrd+=powder.calc_xrd('structures-opt.xyz', CELL, 'XRD-DIFF-OPT-SELF.txt', os.path.abspath('structures-opt-1.cif'), MAX_PROC)
    os.chdir('../')
else:
    topt=0.0
    txrd=0.0

with open('comment.txt','a') as comm:
    comm.write(
        f'''
- calc. name: {job_name};
- system: {project_name};
- default MD time: {MD_steps} steps;
- T = {T} K;
- total of {N_max} iterations;
- greedy exchanges: {exchange_all_pairs}
'''
)


t = time.time()
print(f'\n{"-"*75}\n')
print('Total time: {:.3} h'.format((t-t0)/3600))
print('CP2K/XTB time: {:.3} h'.format(tcalc/3600))
print(f'Optimization: {topt:.3} h')
print(f'Trajectory analysis: {txrd:.3} h')
print('Python overhead: {:.3} s'.format(t-t0-tcalc-(topt+txrd)*3600))

if N_i==N_max:
    print('Normal termination')
else:
    print(f'Error termination: all MDs failed at step {N_i}')