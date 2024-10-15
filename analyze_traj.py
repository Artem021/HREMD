import base_utils
import find_intersections
import os
import multiprocessing as mp
import pandas as pd
from pyxyz import Confpool

MAX_P = 150
STARTING_DATAFILE = '/home/md/md/PES-210_dens/16.07-250_2500-Mn-T600/uniform_polymer_soft.lmps'
# INITIAL_XYZ = '/home/md/md/PES-210_dens/16.07-250_2500-Mn-T600/uniform_polymer.xyz'
# STARTING_XYZ = '/home/md/md/PES-210_dens/16.07-250_2500-Mn-T600/structures.xyz'
STARTING_XYZ = '/home/md/md/PES-210_dens/16.07-250_2500-Mn-T600/alpha_0.7_iter1_restart/traj.xyz'
# STARTING_XYZ = '/home/md/md/PES-185_dens/16.07-250_2500-Mn-T600/test.xyz'
aset = [0.0, 0.05, 0.0625, 0.075, 0.1, 0.3, 0.5, 0.7, 1.0]
aset = [0.7]
p = Confpool()

# proper indexing of dataframe
# for i in range(len(aset)):
    # p.include_from_file(INITIAL_XYZ)
print(f'initial size of confpool: {p.size}')


p.include_from_file(STARTING_XYZ)
print(f'size of confpool: {p.size}')
# 2241
results_final = {i:[] for i in aset}


alphas = {}
for j,alpha in enumerate(aset):
    # alphas[alpha] = [j+i for i in list(range(2250))[::9]][::10] # every 10
    alphas[alpha] = [j+i for i in list(range(2500))][::100] # every 100 (MD 2500)
    # alphas[alpha] = [j+i for i in list(range(2250))[::9]][:11] # first 10
    # alphas[alpha] = [j+i for i in list(range(40))[::9]][::10]
    
# [0] 9 18 27 36 45 54 63 72 81 [90] ...

args=[]
for alp in alphas:
    for i in alphas[alp]:
        argv = ('_default_name.xyz', i, alp, p[i].xyz, p.atom_symbols)
        args.append(argv)
    
print(args)
d = find_intersections.CrossingDetector(datafile=STARTING_DATAFILE,show_progress=False)

# results = [d.analyze_xyz_file(*args[0])]


with mp.Pool(MAX_P) as pool:
    results = pool.starmap(d.analyze_xyz_file, args)
    
for res in sorted(results, key=lambda m: m[1]):
    intersect, indx, alp = res
    
    if intersect is None:
        ninter = float('Nan')
    else:
        ninter = len(intersect[0])
    try:
        results_final[alp].append(ninter)
    except:
        pass
    print(f'      Alpha       Iter      N  ')
    print(f'      {alp}      {indx}   {ninter}')

df = pd.DataFrame.from_dict(results_final)
df.to_excel("PES_210_alpha_0.7_iter1.xlsx")

print(0)








SEQ=False
if SEQ:
    # 1. extract relevant frames from xyz (sorted trajectories from each world)


    # 2. unwrap structures and get number of intersections
    # /usr/bin/python3
    # export LD_LIBRARY_PATH=/home/md/.local/lib/python3.10/site-packages/pyxyz:$LD_LIBRARY_PATH
    # 3. redo energy optimization and trace best structure

    # 4. add energy or N of intersections on graph
    os.chdir(wd)

    # frames = base_utils.readXYZ(fxyz)

    # it=1
    names = ['alpha-0.0.xyz','alpha-0.1.xyz','alpha-0.2.xyz','alpha-0.4.xyz','alpha-0.6.xyz','alpha-0.8.xyz','alpha-1.0.xyz']
    # for j,frame in enumerate(frames):
    #     # if j%len(aset)==0:
    #     #     it+=1
    #     # if it%stride!=0:
    #     #     continue
    #     fname = f'alpha-{aset[j%len(aset)]}.xyz'
    #     names.append(fname)
    #     if os.path.exists(fname):
    #         with open(fname,'a') as fo:
    #             for line in frame:
    #                 fo.write(line)
    #     else:
    #         open(fname,'w').close()


    for file in names:
        d = find_intersections.CrossingDetector(datafile=STARTING_DATAFILE, show_progress=False)
        # with cProfile.Profile() as pr:
        result = d.analyze_xyz_file(file)
        # pr.dump_stats('numba.pstats')

        for i, overlap_data in enumerate(result, start=1):
            if len(overlap_data) == 0:
                print(f"{i}) No ring crossings")
            else:
                print(f"{i}) Ring crossings found: {overlap_data}")
                res = []
                for d in overlap_data:
                    for v in d.values():
                        for elem in v:
                            res.append(elem)
                print(f'Total number of crossings: {len(overlap_data)}')
                print(f'Atom selection for VMD: {" ".join([str(i) for i in set(res)])}')

        r = find_intersections.MoleculeReconstructor(datafile=STARTING_DATAFILE, show_progress=False)
        r.reconstruct_xyz(file, file[:-4]+'-check.xyz')




