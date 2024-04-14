import thread
import multiprocessing as mp
import os, time
LAMMPS_PATH = '/home/artem/LAMMPS/lammps-static/bin/'


alpha_set = ['1.0','0.5','0.0']

worlds = []

alphas = [3.9458734*10**-15, 0.00004, 0.005, 0.01, 0.1, 0.5, 0.9]
alpha = 0.99999999999
for i in range(0,len(alphas)):
    if alpha <= alphas[i]:
        break
    else:
        i = len(alphas)+1
print(i)
alphas.insert(i,alpha)
print(alphas)

exit()
# ---------------- mode 1 ---------------- #

for a in alpha_set:
    w1 = thread.MD_calc(f'/home/artem/LAMMPS_TEST/dbg/alpha-{a}',None,None,[os.path.join(LAMMPS_PATH,'lmp'),'-i','lammps.inp'],
             'lammps.inp','OUT','ERR','traj.xyz','lammps.restart',False,
             float(a),0,'structures.xyz','lammps')
    worlds.append(w1)

worlds[0].print_instance_attributes()

def _run(w):
    w.start_LAMMPS2()
    return time.time(), f'{w.alpha}: {w.last_struc["Energy"]}'

def _set(warr,marr):
    assert len(warr) == len(marr)
    for w in warr:
        pass

# for w in worlds:
#     w.start_LAMMPS2()

for w in worlds:
    print(f'{w.alpha}: {w.last_struc["Energy"]}')
    print(f'length of Eshift: {len(w.ener_shift)} vs MD steps: {w.MD_steps}')

with mp.Pool(3) as pool:
    print(pool.map(_run,worlds))




for w in worlds:
    print(f'{w.alpha}: {w.last_struc["Energy"]}')
    print(f'length of Eshift: {len(w.ener_shift)} vs MD steps: {w.MD_steps}')


