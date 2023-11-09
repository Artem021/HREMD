import os, subprocess, shutil, time
import numpy as np
import xtb_IO
from multiprocessing import Pool
# MAX_PROC = 20
C2_PATH = '/home/artem/'

INPopt = '/home/artem/HREMD/src/cp2k-geo_opt.inp'
TOPopt = '/home/artem/HREMD/src/paracetamol-orthorhombic.prmtop'
TRJopt = 'PAM-pos-1.xyz'

# XYZ = 'structures.xyz'
# WD = '/home/artem/HREMD/multi_opt/'


def run_opt(i,xyz,top=TOPopt,trj=TRJopt,inp=INPopt):
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
                print(f'\t-missing file: {fname}')
                nfail+=1
            ntot+=1

    print(f'\t-Optimization completed: {ntot} of total {nopt} frames successfully optimized!')
    print(f'\t-Number of failed optimizations: {nfail}')
    t=time.time()
    dt=t-t0
    os.chdir('../')
    return dt/60/60



def calc_xrd_matrix():
    pass

def calc_xrd_seq(XYZ, CELL, OUT, REF_CIF,delete_cif=True):
    t0=time.time()
    # 0. convert trajectory to .pdb
    xtb_IO.xyz_to_pdb(XYZ,CELL)
    # 1. create separate .cif files
    if delete_cif:
        xtb_IO.process_traj_xyz(
            fname=XYZ,
            traj_pdb=None,
            dump_pdb=False,
            dump_sep_cif=True,
            cell=CELL,
            ref_cif=REF_CIF
        )
    # 2. calculate difference for each structure# sequential, without full matrix
    CMD_LINE = [os.path.join(C2_PATH,'run_critic.sh')]

    os.chdir('cif')
    cifs = [os.path.abspath(f) for f in os.listdir() if f.endswith('cif')]
    cifs = sorted(cifs, key=lambda x: int(x.split('-')[-1][:-4]))
    cifs.insert(0,os.path.abspath(REF_CIF))
    os.chdir('../')

    diff = []
    istruc=0
    nmax = len(cifs)
    for struc in cifs:
        popen = subprocess.Popen(CMD_LINE, stdout=subprocess.PIPE, stdin=subprocess.PIPE, universal_newlines=True,stderr=subprocess.PIPE,shell=True)
        cmd = f'compare {cifs[0]} {struc}'
        res = popen.communicate(cmd)
        lines = res[0].split('\n')
        # print(lines)
        for line in lines:
            if 'DIFF =' in line and not 'exactly' in line:
                diff.append(line.split()[-1])
                # print(f'{istruc}: {struc}, difference =  {line.split()[-1]}')
                istruc+=1
                # print('Progress: {:.1f}%'.format(istruc/nmax*100),end='\r')
                break
        M = np.array(diff,dtype=float)
    with open(OUT,'w') as fo:
        for i in M:
            fo.write(f'{i:.6f}\n')
    os.chdir('cif')
    shutil.copyfile(f'{XYZ[:-4]}-1.cif',f'../{XYZ[:-4]}-1.cif')
    os.chdir('../')
    t=time.time()
    dt=t-t0
    return dt/60/60
        # return_code = popen.wait()
        # popen.stdout.close()

        # print(return_code)

def calc_difference(struc1,struc2,i=None):
    t0 = time.time()
    proc = subprocess.run(
        os.path.join(C2_PATH,'run_critic.sh'),
        input = f'compare {struc1} {struc2}',
        shell = True,
        text=True,
        capture_output=True,
        )
    t1 = time.time()
    dt=t1-t0
    if proc.returncode!=0:
        return float('nan'),i, dt
    else:
        for line in proc.stdout.split('\n'):
            if '+ DIFF =' in line:
                return float(line.split()[-1]),i,dt

def calc_xrd(XYZ, CELL, OUT, REF_CIF, MAX_PROC):
    # try:
    #     os.mkdir('opt')
    # except:
    #     shutil.rmtree('opt', ignore_errors=False)
    #     os.mkdir('opt')
    # os.chdir('opt')
    # XYZ = 'structures-opt.xyz'
    # # OUT = 'XRD-DIFFERENCE-OPT.txt'
    # OUT = 'XRD-DIFFERENCE-SELF-OPT.txt'
    # # REF_CIF = '/home/artem/HREMD/src/paracetamol-orthorhombic-P1.cif'
    # REF_CIF = '/home/artem/PAM/21.10-12ns/opt/structures-opt-1.cif'
    # CELL={
    #             'a':11.8050,
    #             'b':17.1640,
    #             'c':7.3930,
    #             'alpha':90.00,
    #             'beta':90.00,
    #             'gamma':90.00,
    #             'SpGr': 'Pcab'
    #         }

    # 0. convert trajectory to .pdb
    xtb_IO.xyz_to_pdb(XYZ,CELL)
    # 1. create separate .cif files
    xtb_IO.process_traj_xyz(
        fname=XYZ,
        traj_pdb=None,
        dump_pdb=False,
        dump_sep_cif=True,
        cell=CELL,
        ref_cif=REF_CIF
    )
    os.chdir('cif')
    cifs = [os.path.abspath(f) for f in os.listdir() if f.endswith('cif')]
    cifs = sorted(cifs, key=lambda x: int(x.split('-')[-1][:-4]))
    cifs.insert(0,os.path.abspath(REF_CIF))
    os.chdir('../')

    diff = []
    istruc=0
    nmax = len(cifs)
    times=[]
    indices = []
    t0=time.time()
    with Pool(MAX_PROC) as pool:
        tasks = [(REF_CIF,struc,i) for i,struc in enumerate(cifs)]
        results = [pool.apply_async(calc_difference, t) for t in tasks] #############################
        for r in results:
            res,i,ti = r.get()
            # res = r.get()
            times.append(ti)
            indices.append(i)
            diff.append(res)
            istruc+=1
            # print('Progress: {:.1f}%'.format(istruc/nmax*100),end='\r')

            # print(res)

    # print(diff)
    t1=time.time()
    dt = t1-t0
    # print(f'Total time: {dt/60/60:.2f} h')
    # print(f'Average of 1 pair: {sum(times)/len(times)/60:.2f} min')
    # print(f'Estimated sequential time: {nmax*sum(times)/len(times)/60:.2f} min')
    M = np.array(diff,dtype=float)
    with open(OUT,'w') as fo:
        for i in M:
            fo.write(f'{i:.6f}\n')
    # TODO: make more clear
    os.chdir('cif')
    shutil.copyfile(f'{XYZ[:-4]}-1.cif',f'../{XYZ[:-4]}-1.cif')
    os.chdir('../')
    
    return dt/60/60


# os.chdir('cif')
# shutil.copyfile('structures-1.cif','../structures-1.cif')
# os.chdir('../')







# CELL={
#             'a':11.8050,
#             'b':17.1640,
#             'c':7.3930,
#             'alpha':90.00,
#             'beta':90.00,
#             'gamma':90.00,
#             'SpGr': 'Pcab'
#         }
# os.chdir('/home/artem/PAM/20.10-1ns/')
# calc_xrd_seq('structures.xyz', CELL, OUT = 'XRD-DIFFERENCE.txt', REF_CIF = '/home/artem/HREMD/src/paracetamol-orthorhombic-P1.cif')

# print(0)
# os.chdir('/home/artem/HREMD_algorithm/PAM-crystal/6ns/')
# os.chdir('/home/artem/paracetamol/md/')
# try:
    # os.mkdir('cif')
# except:
    # shutil.rmtree('cif', ignore_errors=True)
    # os.mkdir('cif')
# os.chdir('cif')

# 1. generate separate .cif files
# xtb_IO.process_traj_xyz(
#     # '../PAM-MD-pos-1.xyz',
#     '../structures.xyz',
#     traj_pdb=None,
#     dump_pdb=False,
#     dump_sep_cif=True,
#     cell={
#         'a':11.8050,
#         'b':17.1640,
#         'c':7.3930,
#         'alpha':90.00,
#         'beta':90.00,
#         'gamma':90.00
#     },
#     ref_cif=REF_CIF
# )





# 2. get XRD matrix
# CMD_LINE = [os.path.join(C2_PATH,'critic2')]
# popen = subprocess.Popen(CMD_LINE, stdout=subprocess.PIPE, stdin=subprocess.PIPE, universal_newlines=True,stderr=subprocess.PIPE,shell=True)

# cifs = [os.path.abspath(f) for f in os.listdir() if f.endswith('cif')]
# cifs = sorted(cifs, key=lambda x: int(x.split('-')[-1][:-4]))
# cifs.insert(0,os.path.abspath(REF_CIF))

# cmd = 'compare '+' '.join(cifs)

# res = popen.communicate('compare '+' '.join(cifs))

# lines = res[0].split('\n')
# # print(lines)
# lock=True
# diff = []
# files = []
# for line in lines:
#     if 'DIFF              1               2' in line:
#         Nmax = int(line.split()[-1])
#         lock=False
#         continue
#     if lock:
#         continue
#     if not '.cif' in line:
#         break
#     lspl = line.split()
#     # print(lspl)
#     diff.append(lspl[1:])
#     files.append(lspl[0])
# M = np.array(diff,dtype=float)
# N = np.array(files,dtype=str)

# 3. Write first column to file

# if sequential, M is 1d array





# with open('XRD-2.txt','w') as fo:
#     for i in M:
#         fo.write(f'{i[0]:.6f}\n')

# MT = M.transpose()
# with open('XRD-T.txt','w') as fo:
    # for i in M[0]:
        # fo.write(f'{i:.6f}\n')




