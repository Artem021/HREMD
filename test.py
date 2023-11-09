from multiprocessing import Pool
import os, subprocess, shutil

import xtb_IO

import time, random

INP = '/home/artem/HREMD/src/cp2k-geo_opt.inp'
TOP = '/home/artem/HREMD/src/paracetamol-orthorhombic.prmtop'
TRJ = 'PAM-pos-1.xyz'

# XYZ = 'structures.xyz'
# WD = '/home/artem/HREMD/multi_opt/'


def run_opt(i,xyz,top=TOP,trj=TRJ,inp=INP):
    optdir = f'opt-{i}'
    try:
        os.mkdir(optdir)
    except:
        shutil.rmtree(optdir, ignore_errors=True)
        os.mkdir(optdir)
    shutil.copyfile(top,os.path.join(optdir,os.path.basename(top)))
    shutil.copyfile(inp,os.path.join(optdir,'cp2k.inp'))
    os.chdir(optdir)
    with open('struc.xyz','w') as fo:
        fo.write(xyz)
    # subprocess.call('cp2k.sopt -i cp2k.inp -o cp2k.out', shell=True)
    proc = subprocess.run('cp2k.sopt -i cp2k.inp -o cp2k.out',shell=True)
    if proc.returncode !=0:
        os.chdir('..')
        return f'frame {i}: optimization is down!'
    else:
        with open('struc-opt.xyz','w') as fo:
            xyzopt = xtb_IO.get_frame_xyz(trj,-1)
            fo.write(xyzopt)
    os.chdir('..')
    return f'frame {i}: Normal termination!'

def optimize_all(XYZ,MAX_PROC):
    t0=time.time()
    try:
        os.mkdir('opt')
    except:
        shutil.rmtree('opt', ignore_errors=False)
        os.mkdir('opt')
    os.chdir('opt')
    shutil.copyfile(os.path.join('../',XYZ), XYZ)
    # os.chdir(WD)
    n=0
    structures = []
    with open(XYZ,'r') as traj:
        while True:
            lines = []
            try:
                nat = int(traj.readline())
            except:
                if n==0:
                    print('\t-EOF reached or incorrect number of atoms!')
                else:
                    print(f'\t-Trajectory with {n} structures was read!')
                break
            comm = traj.readline()
            lines.append(f'{nat}\n')
            lines.append(f'Frame {n}\n')
            for i in range(nat):
                lines.append(traj.readline())
            struc = ''.join(lines)
            structures.append(struc)
            n+=1


    print('\t-Start geometry optimization...')

    with Pool(MAX_PROC) as pool:
        tasks = [(i,struc) for i,struc in enumerate(structures)]
        results = [pool.apply_async(run_opt, t) for t in tasks] #############################
        for r in results:
            # print(r.get())
            r.get()



    filenames = [os.path.abspath(f'opt-{i}/struc-opt.xyz') for i in range(len(structures))]

    nopt = 0
    nfail = 0
    ntot = 0
    with open('structures-opt.xyz', 'w') as outfile:
        for fname in filenames:
            try:
                with open(fname) as infile:
                    outfile.write(infile.read()+'\n')
                    nopt+=1
            except:
                print(f'missing file: {fname}')
                nfail+=1
            ntot+=1

    print(f'Optimization completed: {ntot} of total {nopt} frames successfully optimized!')
    print(f'Number of failed optimizations: {nfail}')
    t=time.time()
    dt=t-t0
    return dt/60/60


#! /usr/bin/env python
# import os, subprocess
# import shlex
# import pprint
# # cmd = shlex.split("env -i bash -c 'source /home/artem/cp2k/tools/toolchain/install/setup'")
# cmd = shlex.split("bash -c 'source /home/artem/cp2k/tools/toolchain/install/setup'")
# proc = subprocess.Popen(cmd, stdout = subprocess.PIPE)
# for line in proc.stdout:
#     # print(line)
#     (key, _, value) = line.partition("=")
#     os.environ[key] = value
# #   print(f'{key} ::: {_} ::: {value}')
# proc.communicate()


# pprint.pprint(dict(os.environ))