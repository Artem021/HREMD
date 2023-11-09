import re, os, time
import numpy as np
from multiprocessing import Process
import rmsd
import pyxyz


def combine_xyz(root,logfile,I,Nmax=None):
    N=0
    Nread = 0
    res = []
    fragments = []
    template = 'PAM-pos-1-{n}.xyz'
    with open(logfile,'r') as fi:
        for line in fi:
            res.append([int(i) for i in line.split()])
            N+=1
    print(f'file {logfile} read, {N} H-REMD cycles')
    if not Nmax:
        Nmax=N
    os.chdir(root)
    print(f'enter dir {root}')
    subdirs = [d for d in os.listdir() if os.path.isdir(d)]
    subdirs.sort(key = lambda f: float(f))
    assert(len(subdirs)==len(res[0]))
    while Nread <N-1: # skip last line
        ind = res[Nread][I]
        fname = template.format(n=Nread+1)
        fullname = os.path.join(root,subdirs[ind],fname)
        fragments.append(fullname)
        Nread+=1
    with open(f'replica-{I}.xyz', 'w') as outfile:
        for fname in fragments:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)
    return fragments

# res = combine_xyz('/home/artem/HREMD_algorithm/cp2k/ATROPOISOMERS-VEL-SOFT-8',
#                   '/home/artem/HREMD_algorithm/cp2k/ATROPOISOMERS-VEL-SOFT-8/exchanges.log',
#                   0)
# print(res)

# https://stackoverflow.com/questions/20305272/dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python

def get_dihedral2(coord):
    p = np.array(coord,dtype=float)
    b = p[:-1] - p[1:]
    b[0] *= -1
    v = np.array( [ v - (v.dot(b[1])/b[1].dot(b[1])) * b[1] for v in [b[0], b[2]] ] )
    # Normalize vectors
    v /= np.sqrt(np.einsum('...i,...i', v, v)).reshape(-1,1)
    b1 = b[1] / np.linalg.norm(b[1])
    x = np.dot(v[0], v[1])
    m = np.cross(v[0], b1)
    y = np.dot(m, v[1])
    return np.degrees(np.arctan2( y, x ))


def get_dihedral(coord):
    # coord = [[x,y,z],[x,y,z],[x,y,z],[x,y,z]]
    m1,m2,m3,m4 = coord
    a1 = np.array(m1,dtype=float)
    a2 = np.array(m2,dtype=float)
    a3 = np.array(m3,dtype=float)
    a4 = np.array(m4,dtype=float)
    v1 = a3 - a1
    v2 = a2 - a1
    cp1 = np.cross(v1, v2)
    v3 = a4 - a2
    v4 = a3 - a2
    cp2 = np.cross(v3, v4)
    cosine_angle = np.dot(cp1, cp2) / (np.linalg.norm(cp1) * np.linalg.norm(cp2))
    # angle = np.arccos(cosine_angle)
    angle = np.arccos(np.clip(cosine_angle, -1, 1))
    return np.degrees(angle).round(3)

def check_xyz(xyzname,res = [],dh = [16,15,14,50],filter=False):
    if filter:
        rs = 0
        xyzname2 = xyzname[:-4]+'-filter'+xyzname[-4:]
        fo = open(xyzname2,'w')
        open('group1.xyz','w').close()
        open('group2.xyz','w').close()
        open('group3.xyz','w').close()
    with open(xyzname,'r') as fi:
        while True:
            xyzdh = [0,0,0,0]
            xyzlines = []
            try:
                line1 = fi.readline()
                xyzlines.append(line1)
                nat = int(line1)
            except: # EOF reached
                break
            comm = fi.readline()
            xyzlines.append(comm)
            for i in range(1,nat+1):
                line = fi.readline()
                xyzlines.append(line)
                if i in dh:
                    xyzdh[dh.index(i)] = line.split()[1:]
            # dhval = get_dihedral(xyzdh)
            dhval = get_dihedral2(xyzdh)
            res.append(dhval)
            if filter:
                if 0<= dhval <= 80:
                    with open('group1.xyz','a') as g1:
                        g1.write(''.join(xyzlines))
                elif 80<= dhval <= 135:
                    with open('group2.xyz','a') as g2:
                        g2.write(''.join(xyzlines))
                elif 135<= dhval <= 180:
                    with open('group3.xyz','a') as g3:
                        g3.write(''.join(xyzlines))
                if dhval > 0:
                    fo.write(''.join(xyzlines))
                else:
                    rs+=1
                    print('struc removed')
    if filter:
        print(f'Total of {rs} removed structures')
        fo.close()
    return res

def rmsd_xyz(xyzname,refcrd,res = []):
    with open(xyzname,'r') as fi:
        while True:
            xyz = []
            try:
                nat = int(fi.readline())
            except: # EOF reached
                break
            comm = fi.readline()
            for i in range(nat):
                line = fi.readline()
                xyz.append(line.split()[1:])
            crd = np.array(xyz,dtype=float)
            val = rmsd.kabsch_rmsd(refcrd, crd, translate=True)
            res.append(val)
    return res
            
def rmsd_pyxyz(fnames,refname):
    rmsd_values =[]
    p = pyxyz.Confpool()
    p.include_from_file(refname)
    print(f'Reading input files:')
    if isinstance(fnames, str):
        p.include_from_file(fnames)
        print(f'File read successfully')
    else:
        n=0
        nmax = len(fnames)
        for f in fnames:
            p.include_from_file(f)
            n+=1
            print('Progress: {:.1f}%'.format(n/nmax*100),end='\r')
        print(f'Total number of .xyz files read: {nmax}')
    # p.generate_connectivity(0, mult=1.3, ignore_elements=['HCarbon'], sdf_name='connectivity_check.sdf')
    p.generate_connectivity(0, mult=1.3, ignore_elements=['HCarbon'])
    niso = p.generate_isomorphisms()
    ref = p[0]
    n=0
    nmax = len(p)-1
    for i,struc in enumerate(p):
        val, U, center = struc.rmsd(ref)
        if i==0: # skip 0th structure
            continue
        rmsd_values.append(val)
        n+=1
        print('Progress: {:.1f}%'.format(n/nmax*100),end='\r')
    print('Success!')
    return rmsd_values

def check_dir(dirname, refname=None, refcrd=None, log='results.txt', do_rmsd = True, write=False, use_pyxyz=True):
    alpha = os.path.basename(dirname)
    files = [os.path.join(dirname,f) for f in os.listdir(dirname) if f.endswith('xyz')]
    files.sort(key=lambda f: os.path.getmtime(f))
    files.pop() # last modified .xyz file is always input-structure
    nf = len(files)
    n = 0
    res = []
    for f in files:
        if do_rmsd:
            if use_pyxyz:
                res = rmsd_pyxyz(files,refname)
                return (alpha,res)
            assert refcrd!=None
            res = rmsd_xyz(f,refcrd,res)
            # res = rmsd_pyxyz(f,refname,res)
        else:    
            res = check_xyz(f,res)
        n+=1
        print('Progress: {:.1f}%'.format(n/nf*100),end='\r')
    if write:
        with open(log,'w') as fo:
            for val in res:
                fo.write(f'{val}\n')
    else:
        return (alpha,res)

def calc_rmsd(refxyz,maindir,logname,do_rmsd=True):
    _, refcrd = rmsd.get_coordinates(refxyz, "xyz")
    t0 = time.time()
    os.chdir(maindir)
    resdirs = [d for d in os.listdir(maindir) if os.path.isdir(os.path.join(maindir,d))]
    resdirs.sort(key=lambda f: float(f))
    finres = []
    idir = 0
    ndir = len(resdirs)
    for d in resdirs:
        # print(f'Processing directory: {d}',end='\r')
        print(f'Processing directory: {d}')
        alpha,res = check_dir(os.path.join(maindir,d),refcrd=refcrd,refname=refxyz,do_rmsd=do_rmsd)
        finres.append(res)
    n = max([len(a) for a in finres])
    with open(logname,'w') as fo:
        fo.write(''.join([alpha.center(25) for alpha in resdirs])+'\n') # header
        for i in range(n):
            fo.write(''.join([f'{a[i]}'.center(25) for a in finres])+'\n')
    t = time.time()
    print(f'Total time: {(t-t0)/60} m')

def get_structure_xyz(xyz):
    _, crd = rmsd.get_coordinates(xyz, "xyz")
    return crd

if __name__ == '__main__':

    os.chdir('/home/artem/HREMD_algorithm/cp2k/rmsd/crest_results')
    # val1 = rmsd_pyxyz(f'crest_conformers-1.xyz','first.xyz')
    # val2 = rmsd_pyxyz(f'crest_conformers-2.xyz','first.xyz')
    val1 = check_xyz(f'crest_conformers-1.xyz',res=[])
    val2 = check_xyz(f'crest_conformers-2.xyz',res=[])
    val3 = check_xyz(f'/home/artem/HREMD_algorithm/cp2k/rmsd/crest_results/structures-2-filter.xyz',res=[])
    # val4 = check_xyz(f'/home/artem/HREMD_algorithm/cp2k/rmsd/opt/vel_soft-8/structures-opt.xyz',res=[])
    # val4 = check_xyz(f'/home/artem/HREMD_algorithm/cp2k/rmsd/opt/vel_soft-8-xtbopt/structures-opt-gfnff.xyz',res=[])
    val4 = check_xyz(f'/home/artem/HREMD_algorithm/cp2k/rmsd/crest_results/structures-2.xyz',res=[])
    val5 = check_xyz(f'/home/artem/HREMD_algorithm/cp2k/rmsd/opt/vel_soft-8-xtbopt/structures-opt-gfnff-filter.xyz',
                    res=[],
                    #  dh=[49,24,33,34],
                    filter=True)
    with open('dihedral-hremd-2-filter.txt','w') as fo:
        for i in val3:
            fo.write(f'{i}\n'.ljust(15))

    with open('improper-hremd-1.txt','w') as fo:
        for i in val5:
            fo.write(f'{i}\n'.ljust(15))
    with open('dihedral-hremd-2.txt','w') as fo:
        for i in val4:
            fo.write(f'{i}\n'.ljust(15))

    with open('dihedral-crest-1.txt','w') as fo:
        for i in val1:
            fo.write(f'{i}\n'.ljust(15))
    with open('dihedral-crest-2.txt','w') as fo:
        for i in val2:
            fo.write(f'{i}\n'.ljust(15))
    print(val1)
    #
    #print(f'RMSD of candidate and isomer 2: {val1}')
    # print(f'RMSD of candidate and isomer 2: {val2}')
    # val3 = rmsd_pyxyz(f'atropo1-{prog}opt.xyz',f'atropo2-{prog}opt.xyz')
    # print(f'RMSD of isomers 1 and 2: {val3}')

    # val4 = rmsd_pyxyz('/home/artem/HREMD_algorithm/cp2k/rmsd/opt/atropo-1/atropo-1-opt-pos-1.xyz',f'atropo1-{prog}opt.xyz')
    # print(f'RMSD of opt. trajectory of isomer 1: {val4}')

    # val5 = rmsd_pyxyz('/home/artem/HREMD_algorithm/cp2k/rmsd/opt/atropo-2/atropo-2-opt-pos-1.xyz',f'atropo2-{prog}opt.xyz')
    # print(f'RMSD of opt. trajectory of isomer 2: {val5}')

    # val6 = rmsd_pyxyz('/home/artem/HREMD_algorithm/cp2k/rmsd/opt/candidate-1/cand-1-opt-pos-1.xyz',f'cand1-{prog}opt.xyz')
    # print(f'RMSD of opt. trajectory of candidate: {val6}')

    # val7 = rmsd_pyxyz('/home/artem/HREMD_algorithm/cp2k/rmsd/relax/atropo-2-relax-pos-1.xyz','second.xyz')
    # print(f'RMSD of relax. trajectory of isomer 2: {val7}')

    # # exit()
    # dirname = '/home/artem/HREMD_algorithm/cp2k/ATROPOISOMERS-SHAKE/1.000e+00'
    # alpha = os.path.basename(dirname)
    # files1 = [os.path.join(dirname,f) for f in os.listdir(dirname) if f.endswith('xyz')]
    # files1.sort(key=lambda f: os.path.getmtime(f))
    # files1.pop() # last modified .xyz file is always input-structure

    # dirname = '/home/artem/HREMD_algorithm/cp2k/ATROPOISOMERS-SHAKE/0.000e+00'
    # alpha = os.path.basename(dirname)
    # files2 = [os.path.join(dirname,f) for f in os.listdir(dirname) if f.endswith('xyz')]
    # files2.sort(key=lambda f: os.path.getmtime(f))
    # files2.pop() # last modified .xyz file is always input-structure


    # calc_rmsd('/home/artem/HREMD_algorithm/cp2k/rmsd/second.xyz','/home/artem/HREMD_algorithm/cp2k/ATROPOISOMERS-SHAKE/','rmsd-alt.txt')
    # calc_rmsd('/home/artem/HREMD_algorithm/cp2k/rmsd/first.xyz','/home/artem/HREMD_algorithm/cp2k/ATROPOISOMERS-SHAKE/','dihedral.txt',do_rmsd=False)
    calc_rmsd('/home/artem/HREMD_algorithm/cp2k/rmsd/atropo2-cp2kopt.xyz','/home/artem/HREMD_algorithm/cp2k/ATROPOISOMERS-NO-SHAKE/','rmsd-fix.txt')
    calc_rmsd('/home/artem/HREMD_algorithm/cp2k/rmsd/atropo2-cp2kopt.xyz','/home/artem/HREMD_algorithm/cp2k/ATROPOISOMERS-SHAKE/','rmsd-fix.txt')
    # alpha, res = check_dir('/home/artem/HREMD_algorithm/cp2k/ATROPOISOMERS-SHAKE/0.000e+00',
    #                        refname='/home/artem/HREMD_algorithm/cp2k/rmsd/second.xyz',
    #                        do_rmsd = False)
    # print(alpha)