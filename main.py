import os,math,random,shutil,time,queue
from multiprocessing import Process,Queue, Pool
import thread
import results as postprocess
import powder

note = '''
- T-REMD;
- default MD time: 1000 steps;
- FIXED BUGS:
    1. when smooth alpha affected all replicas
    2. in `update_shared_xyz`: struc_before_exchange --> struc_after_exchange
- compound: Paracetamol;
- periodic calculation in NVT with orthorhombic cell;
- T = {1450 : ... : 298} K;
- start from monoclinic structure optimized in orthorhombic cell (guess-1.xyz);
- total of 500*6 iterations (6 ns);
- greedy exchanges: YES;
- smooth alpha: YES,
    alpha[max] = 0.05,
    rise if alpha <0.05;
- external potential: YES,
    everywhere except a>=0.1,
    formula: C*((DIM)^2)^4, where C=10**-14; DIM = X,Y,Z,
    potential affects all 160 atoms;
- 'shaking' cells: NO,
    dr = 20,
    C=10**-10;
- velocity softening: YES,
    100 steps,
    default alpha and delta;
- default 1-4 interactions (1.0 electrostatic, 1.0 VdW);
- full contributions from torsions and improper torsions
'''

#----------project general settings----------#

JOB_NAME = '08.11-6ns-TREMD'
TOP = 'paracetamol-orthorhombic.prmtop'
PROJECT_NAME = 'PAM'
XYZ = 'guess-1.xyz'
if XYZ[-3:]=='pdb':
    XYZ_EXT = 'PDB'
elif XYZ[-3:]=='xyz':
    XYZ_EXT = 'XYZ'
else:
    print(f'unknown format of coordinate file: {XYZ}')
    raise ValueError

#----------MD settings----------#

alpha_vs = 0.15 # default: 0.15
delta_vs = 0.10 # default: 0.10
steps_vs = 100 # default: 40
md_steps = 1000

TEMPERATURE = 298 # default: 298

PERIODIC = True

CELL={
            'a':11.8050,
            'b':17.1640,
            'c':7.3930,
            'alpha':90.00,
            'beta':90.00,
            'gamma':90.00,
            'SpGr': 'Pcab'
        }


#----------HREMD settings----------#

# alpha_set = [0.0, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
alpha_set = [1450, 1400, 1300, 1200, 1150, 1100, 1050, 1000, 900, 800, 700, 600, 500, 400, 298]
# alpha_set = [700, 650, 600, 575, 550, 525, 500, 475, 450, 425, 400, 375, 350, 325, 298]
T_REMD = True
if T_REMD:
    print(f'{"="*50}\n{" Do classic T-REMD! ".center(50,"=")}\n{"="*50}')

NMAX = 500*6 # Total number of H-REMD cycles (500 for total MD time of 1 ns)

MAX_PROC = 22

BETA = 1059.5151978
kB = 3.167*10**-6 # Boltzmann k in Ha/K

XYZ_SHARED = 'structures.xyz'

EXCHANGE_STAT = 'exchanges.log'

ENERGY_STAT = 'energy.log'

RMSD_REF = None

GREEDY_EXCHANGE = True

USE_LAST_STRUC = False # if lowest energy structure is in first half of trajectory, the last structure selected

EXTERNAL_POT = True

ATOMS_EXT = "1..160" # atoms affected by external potential

USE_RMSD = False # TODO

CP2K_SHAKE_CELL = False

PARM = {
    "GLOBAL": [
        {
            "PRINT_LEVEL": [
                "MEDIUM"
            ],
            "PROJECT": [
            f"{PROJECT_NAME}"
            ],
            "RUN_TYPE": [
                "MD"
            ],
            "SEED": [
                "1642"
            ]
        }
    ],
    "FORCE_EVAL": [
        {
            "EXTERNAL_POTENTIAL": [
                {
                    "ATOMS_LIST": [
                        f"{ATOMS_EXT}"
                    ],
                    "FUNCTION": [
                        "0.001*0.00000000001*(Z^2)^4"
                    ]
                },
                {
                    "ATOMS_LIST": [
                        f"{ATOMS_EXT}"
                    ],
                    "FUNCTION": [
                        "0.001*0.00000000001*(X^2)^4"
                    ]
                },
                {
                    "ATOMS_LIST": [
                        f"{ATOMS_EXT}"
                    ],
                    "FUNCTION": [
                        "0.001*0.00000000001*(Y^2)^4"
                    ]
                }
            ],
            "METHOD": [
                "FIST"
            ],
            "MM": [
                {
                    "FORCEFIELD": [
                        {
                            # "IMPROPER": [
                            #     {
                            #         "ATOMS": [
                            #             "27 25 26 28"
                            #         ],
                            #         "K": [
                            #             "100"
                            #         ],
                            #         "PHI0": [
                            #             "0"
                            #         ]
                            #     },
                            #     {
                            #         "ATOMS": [
                            #             "49 24 33 34"
                            #         ],
                            #         "K": [
                            #             "100"
                            #         ],
                            #         "PHI0": [
                            #             "0"
                            #         ]
                            #     },
                            #     {
                            #         "ATOMS": [
                            #             "2 10 3 4"
                            #         ],
                            #         "K": [
                            #             "100"
                            #         ],
                            #         "PHI0": [
                            #             "0"
                            #         ]
                            #     }
                            # ],
                            "DO_NONBONDED": [
                                "T"
                            ],
                            "ALL_EI_SCALE": [
                                "1.0"
                            ],
                            "ALL_VDW_SCALE": [
                                "1.0"
                            ],
                            "EI_SCALE14": [
                                "1.0"
                            ],
                            "VDW_SCALE14": [
                                "1.0"
                            ],
                            "FORCE_SCALE": [
                                "1.0"
                            ],
                            "PARMTYPE": [
                                "AMBER"
                            ],
                            "PARM_FILE_NAME": [
                                f"{TOP}"
                            ],
                            "SHIFT_CUTOFF": [
                                "F"
                            ],
                            "SPLINE": [
                                {
                                    "EMAX_SPLINE": [
                                        "1.0E+08"
                                    ],
                                    "RCUT_NB": [
                                        "[angstrom] 1.0E+01"
                                    ]
                                }
                            ]
                        }
                    ],
                    "POISSON": [
                        {
                            "EWALD": [
                                {
                                    "ALPHA": [
                                        ".35" if PERIODIC else ".40"
                                    ],
                                    "EWALD_TYPE": [
                                        "EWALD" if PERIODIC else "NONE"
                                    ],
                                    "GMAX": [
                                        "25" if PERIODIC else "80"
                                    ]
                                }
                            ],
                            "PERIODIC": [
                                "XYZ" if PERIODIC else "NONE"
                            ]
                        }
                    ]
                }
            ],
            "RESCALE_FORCES": [
                {
                    "MAX_FORCE": [
                        "10"
                    ]
                }
            ],
            "STRESS_TENSOR": [
                "ANALYTICAL" if PERIODIC else "NONE"
            ],
            "SUBSYS": [
                {
                    "CELL": [
                        {
                            "ABC": [
                                "[angstrom] " + f"{CELL['a']:.4f}  {CELL['b']:.4f}  {CELL['c']:.4f}" if PERIODIC else "100 100 100"
                            ],
                            "ALPHA_BETA_GAMMA": [
                                f"{CELL['alpha']:.2f} {CELL['beta']:.2f} {CELL['gamma']:.2f}" if PERIODIC else "90 90 90"
                            ],
                            "PERIODIC": [
                                "XYZ" if PERIODIC else "NONE"
                            ]
                        }
                    ],
                    "TOPOLOGY": [
                        {
                            "CENTER_COORDINATES": [
                                {
                                    "CENTER_POINT": [
                                        "0 0 0"
                                    ]
                                }
                            ],
                            "CONN_FILE_FORMAT": [
                                "AMBER"
                            ],
                            "CONN_FILE_NAME": [
                                f"{TOP}"
                            ],
                            "COORD_FILE_FORMAT": [
                                "XYZ" #TODO: add list of known file types
                            ],
                            "COORD_FILE_NAME": [
                                f"{XYZ}"
                            ]
                        }
                    ]
                }
            ]
        }
    ],
    "MOTION": [
        {
            "MD": [
                {
                    "VELOCITY_SOFTENING": [
                        {
                            "STEPS": [
                                f"{steps_vs}"
                            ],
                            "ALPHA": [
                                f"{alpha_vs}"
                            ],
                            "DELTA ": [
                                f"{delta_vs}"
                            ]
                        }
                    ],
                    "ENSEMBLE": [
                        "NVT"
                    ],
                    "STEPS": [
                        f"{md_steps}"
                    ],
                    "TEMPERATURE": [
                        f"{TEMPERATURE}"
                    ],
                    "THERMOSTAT": [
                        {
                            "CSVR": [
                                {
                                    "TIMECON": [
                                        "[fs] 0.5"
                                    ]
                                }
                            ],
                            "REGION": [
                                "GLOBAL"
                            ],
                            "TYPE": [
                                "CSVR"
                            ]
                        }
                    ],
                    "TIMESTEP": [
                        "[fs] 2.0"
                    ]
                }
            ],
            "PRINT": [
                {
                    "RESTART": [
                        {
                            "BACKUP_COPIES": [
                                "0"
                            ],
                            "EACH": [
                                {
                                    "MD": [
                                        "0"
                                    ]
                                }
                            ]
                        }
                    ],
                    "RESTART_HISTORY": [
                        {
                            "EACH": [
                                {
                                    "MD": [
                                        "0"
                                    ]
                                }
                            ]
                        }
                    ],
                    "TRAJECTORY": [
                        {
                            "EACH": [
                                {
                                    "MD": [
                                        "5"
                                    ]
                                }
                            ],
                            "FORMAT": [
                                "XYZ" #TODO: add PDB parser
                            ]
                        }
                    ]
                }
            ]
        }
    ]
}


#----------set directory----------#

Nxyz = 0
wd = '../' + PROJECT_NAME

try:
    os.mkdir(wd)
except:
    pass
shutil.copyfile(os.path.join('src', XYZ), os.path.join(wd, XYZ))
shutil.copyfile(os.path.join('src', TOP), os.path.join(wd, TOP))
os.chdir(wd)
try:
    os.mkdir(JOB_NAME)
except:
    shutil.rmtree(JOB_NAME, ignore_errors=True)
    os.mkdir(JOB_NAME)
try:
    shutil.copyfile(XYZ_SHARED,os.path.join(JOB_NAME,XYZ_SHARED))
    with open(XYZ_SHARED,'r') as fi: # get N of initial structures
        natoms = int(fi.readline())
        nlines = sum(1 for i in fi if i != '\n') + 1
        Nxyz = int(nlines/(2+natoms))
    print(f'{Nxyz} initial structures were read')
except:
    Nxyz = 0
    print('initial structures were not provided')
shutil.copyfile(XYZ,os.path.join(JOB_NAME,XYZ))
shutil.copyfile(TOP,os.path.join(JOB_NAME,TOP))
os.chdir(JOB_NAME)

#----------create Worlds----------#
worlds = []
for alp in alpha_set:
    world = thread.World(
        BASE_DIR = os.path.abspath(os.getcwd()),
        alpha = alp,
        T_REMD = T_REMD,
        TOP = TOP,
        XYZ = XYZ,
        # XYZ0 = XYZ,
        XYZ_FOUND=XYZ_SHARED,
        ARGS = PARM, # __init__ takes deepcopy of dict
        USE_LAST_STRUC = USE_LAST_STRUC
    )
    worlds.append(world)
    if world.alpha < 0.001 and CP2K_SHAKE_CELL: # TODO: make it more clear
        world.shake = True
    if not EXTERNAL_POT or world.alpha >= 0.1:
        world.remove_ext_pot()

# start of loop
exchange_set = ['{:.3e}'.format(i) for i in alpha_set]
exchange_ord = {alp:exchange_set.index(alp) for alp in exchange_set}
with open(EXCHANGE_STAT,'w') as exc_stat:
    exc_stat.write('   '.join([str(exchange_ord[a]) for a in exchange_set]) + '\n')
with open(ENERGY_STAT,'w') as e_stat:
    e_stat.write('Type'.center(25) + ''.join([a.center(25) for a in exchange_set]) + '\n')
with open('comments.txt','w') as comm:
    comm.write(note)

world1 = worlds[-1]
# world1.remove_ext_pot()
N = 0
txtb = 0
t0 = time.time()
# if CP2K_SHAKE_CELL:
#     for w in worlds:
#         w.set_pot()

while N < NMAX:
    # 1. parallel launch of MDs
    procs = []
    main_queue = Queue()
    results = []
    for i,w in enumerate(worlds):
        if w.shake:
            w.shake_cell(N)
            # w.set_pot()
        calc = thread.MD_calc(w.XYZ,w.ARGS,w.CMD_LINE,w.INP,w.OUT,w.ERR,w.TRJ,w.alpha,Nxyz,w.XYZ_FOUND)
        proc = Process(target = calc.start_XTB,args=(main_queue,i))
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
    txtb+=max([r['time']for r in results])
    # if all MDs failed, we exit
    if all([w.error for w in worlds]):
        print(f'All calculations are down, exiting (cycle {N})')
        break
    # 2. Metropolis–Hastings part
    print(f'Cycle {N}')
    for w in worlds:
        w.exchanged = False
    for i in range(len(worlds)-1):
        w1 = worlds[i]
        w2 = worlds[i+1]
        if (w1.error | w2.error):
            continue
        if w1.exchanged: # current implementation permits non-neighbors to exchange
            if GREEDY_EXCHANGE:
                struc1 = w1.get_structure_after()  
            else:
                continue
        else:
            struc1 = w1.get_structure_before()
        if w2.exchanged:
            if GREEDY_EXCHANGE:
                struc2 = w2.get_structure_after()  
            else:
                continue
        else:
            struc2 = w2.get_structure_before()
        alp1 = w1.alpha
        alp2 = w2.alpha
        e1 = struc1['Energy']
        e2 = struc2['Energy']
        if T_REMD:
            beta1 = (kB*alp1)**-1
            beta2 = (kB*alp2)**-1
            rel_e = -e1*beta1 -e2*beta2 + e2*beta1 + e1*beta2 # delta = (e1-e2)*(alp2-alp1)
            delta = max(-400,min(400,rel_e)) # returns -400 if rel_e==nan
        else:
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
        if USE_RMSD:
            struc1 = w1.get_structure_before()
            struc2 = w2.get_structure_before()
            val1 = struc1['RMSD']
            val2 = struc2['RMSD']
            exchange = (val1 < val2)
            print('Using RMSD criteria instead..')
            print(f'RMSD of {alp1}: {val1} | RMSD of {alp2}: {val2}')
        if exchange:
            w1.set_structure_after(struc2)
            w2.set_structure_after(struc1)
            w1.update_xyz()
            w2.update_xyz()
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
            w1.update_xyz()
            w2.update_xyz()
        # print('=AFTER EXCHANGE ATTEMPT=')
        # print(f'struc1, before: {w1.get_structure_before()}')
        # print(f'struc2, before: {w2.get_structure_before()}')
        # print(f'struc1, after: {w1.get_structure_after()}')
        # print(f'struc2, after: {w2.get_structure_after()}')




    with open(EXCHANGE_STAT,'a') as exc_stat:
        exc_stat.write('   '.join([str(exchange_ord[a]) for a in exchange_set]) + '\n')
    with open(ENERGY_STAT,'a') as e_stat:
        e_stat.write('ENERGY'.center(25) + ''.join(['{:.8e}'.format(w.struc_before_exchange['Energy']).center(25) for w in worlds]) + '\n')
        e_stat.write('VBIAS'.center(25) + ''.join(['{:.8e}'.format(w.struc_before_exchange['Vbias']).center(25) for w in worlds]) + '\n')
        if USE_RMSD:
            e_stat.write('RMSD'.center(25) + ''.join(['{:.8e}'.format(w.struc_before_exchange['RMSD']).center(25) for w in worlds]) + '\n')
    # 3. update shared XYZ
    world1.update_shared_xyz()
    Nxyz+=1
    
    # 4. prepare for next step
    for w in worlds:
        if w.error:
            #TODO: update .xyz with previos/from another world
            continue
        w.rename_trj_file()
    N+=1

    # 5. change order of exchange attempts
    # worlds.reverse() # TODO: better understand this feature

# postprocessing
# if RMSD_REF != None:
#     postprocess.calc_rmsd(RMSD_REF,os.path.abspath(os.getcwd()),'rmsd.txt')
#     postprocess.calc_rmsd(RMSD_REF,os.path.abspath(os.getcwd()),'dihedral.txt',do_rmsd=False)




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


t = time.time()
print('Total time: {:.3} h'.format((t-t0)/3600))
print('CP2K/XTB time: {:.3} h'.format(txtb/3600))
print(f'Optimization: {topt:.3} h')
print(f'Trajectory analysis: {txrd:.3} h')
print('Python overhead: {:.3} s'.format(t-t0-txtb-(topt+txrd)*3600))

if N==NMAX:
    print('Normal termination')
else:
    print(f'Error termination: all MDs failed at step {N}')