import os, re, collections, shutil
import numpy as np

TEST = False

DEFAULT_PARM_CP2K = {
    'GLOBAL': {
        'SEED': '1642',
        'PROJECT': 'default', 
        'PRINT_LEVEL': 'MEDIUM', 
        'RUN_TYPE': 'MD'
    },
    'FORCE_EVAL': {
        'RESCALE_FORCES': {
            'MAX_FORCE': '10'
        }, 
        'STRESS_TENSOR': 'ANALYTICAL', 
        'METHOD': 'FIST', 
        'MM': {
            'FORCEFIELD': {
                'FORCE_SCALE': '1.0', 
                'EI_SCALE14': '1.0', 
                'ALL_VDW_SCALE': '1.0', 
                'ALL_EI_SCALE': '1.0', 
                'SHIFT_CUTOFF': 'F', 
                'PARMTYPE': 'AMBER', 
                'PARM_FILE_NAME': 'default.prmtop', 
                'SPLINE': {
                    'EMAX_SPLINE': '1.0E+08', 
                    'RCUT_NB': '[angstrom] 1.0E+01'
                }
            }, 
            'POISSON': {
                'EWALD': {
                    'EWALD_TYPE': 'NONE', 
                    'ALPHA': '.40', 
                    'GMAX': '80'
                }, 
                'PERIODIC': 'NONE'
            }
        }, 
        'SUBSYS': {
            'CELL': {
                'ABC': '[angstrom] 100 100 100', 
                'ALPHA_BETA_GAMMA': '90 90 90', 
                'PERIODIC': 'NONE'
                }, 
            'TOPOLOGY': {
                'CENTER_COORDINATES': {}, 
                'CONN_FILE_FORMAT': 'AMBER', 
                'CONN_FILE_NAME': 'default.prmtop', 
                'COORD_FILE_FORMAT': 'XYZ', 
                'COORD_FILE_NAME': 'default.XYZ'
            }
        }
    },
    'MOTION': {
        'MD': {
            'ENSEMBLE': 'NVT', 
            'TIMESTEP': '[fs] 2.0', 
            'STEPS': '1000', 
            'TEMPERATURE': '298', 
            'THERMOSTAT': {
                'REGION': 'GLOBAL', 
                'TYPE': 'CSVR', 
                'CSVR': {
                    'TIMECON': '[fs] 0.5'
                }
            }
        }, 
        'PRINT': {
            'RESTART': {
                'BACKUP_COPIES': '0', 
                'EACH': {
                    'MD': '0'
                }
            }, 
            'TRAJECTORY': {
                'FORMAT': 'XYZ', 
                'EACH': {
                    'MD': '5'
                }
            }, 
            'RESTART_HISTORY': {
                'EACH': {
                    'MD': '0'
                }
            }
        }
    }
}

DEFAULT_PARM = {
    '$md':{'restart':'false','hmass':4,'time':10.00,'temp':300.00,'step':2.00,'shake':0,'dump':10.00},
    '$metadyn':{'save':5,'kpush':0.042000,'alp':1.300000,'alpha':0.0},
    '$constrain':{'force constant':0.25,'all bonds':'T','reference':'reference.xyz'},
    '$wall':{'potential':'logfermi','sphere':'auto,all'}
}

SPECIAL_SIGNS = {'elements':':','distance':':','sphere':':'}

def get_fractional_to_cartesian_matrix(a, b, c, alpha, beta, gamma,
                                       angle_in_degrees=True):
    """
    Return the transformation matrix that converts fractional coordinates to
    cartesian coordinates.

    Parameters
    ----------
    a, b, c : float
        The lengths of the edges.
    alpha, gamma, beta : float
        The angles between the sides.
    angle_in_degrees : bool
        True if alpha, beta and gamma are expressed in degrees.

    Returns
    -------
    r : array_like
        The 3x3 rotation matrix. ``V_cart = np.dot(r, V_frac)``.

    """
    if angle_in_degrees:
        alpha = np.deg2rad(alpha)
        beta = np.deg2rad(beta)
        gamma = np.deg2rad(gamma)
    cosa = np.cos(alpha)
    sina = np.sin(alpha)
    cosb = np.cos(beta)
    sinb = np.sin(beta)
    cosg = np.cos(gamma)
    sing = np.sin(gamma)
    volume = 1.0 - cosa**2.0 - cosb**2.0 - cosg**2.0 + 2.0 * cosa * cosb * cosg
    volume = np.sqrt(volume)
    r = np.zeros((3, 3))
    r[0, 0] = a
    r[0, 1] = b * cosg
    r[0, 2] = c * cosb
    r[1, 1] = b * sing
    r[1, 2] = c * (cosa - cosb * cosg) / sing
    r[2, 2] = c * volume / sing
    return r


def get_cartesian_to_fractional_matrix(a, b, c, alpha, beta, gamma,
                                       angle_in_degrees=True):
    """
    Return the transformation matrix that converts cartesian coordinates to
    fractional coordinates.

    Parameters
    ----------
    a, b, c : float
        The lengths of the edges.
    alpha, gamma, beta : float
        The angles between the sides.
    angle_in_degrees : bool
        True if alpha, beta and gamma are expressed in degrees.

    Returns
    -------
    r : array_like
        The 3x3 rotation matrix. ``V_frac = np.dot(r, V_cart)``.

    """
    if angle_in_degrees:
        alpha = np.deg2rad(alpha)
        beta = np.deg2rad(beta)
        gamma = np.deg2rad(gamma)
    cosa = np.cos(alpha)
    sina = np.sin(alpha)
    cosb = np.cos(beta)
    sinb = np.sin(beta)
    cosg = np.cos(gamma)
    sing = np.sin(gamma)
    volume = 1.0 - cosa**2.0 - cosb**2.0 - cosg**2.0 + 2.0 * cosa * cosb * cosg
    volume = np.sqrt(volume)
    r = np.zeros((3, 3))
    r[0, 0] = 1.0 / a
    r[0, 1] = -cosg / (a * sing)
    r[0, 2] = (cosa * cosg - cosb) / (a * volume * sing)
    r[1, 1] = 1.0 / (b * sing)
    r[1, 2] = (cosb * cosg - cosa) / (b * volume * sing)
    r[2, 2] = sing / (c * volume)
    return r



def read_cp2k(fname,fi=None):
    sect = {}
    if fi == None:
        fi = open(fname,'r')
    while True:
        line = fi.readline()
        if line == '\n':
            continue
        if 'END' in line or line == '':
            return sect
        line = line.strip()
        if '&' in line:
            sctname = line[1:]
            subsect = read_cp2k(fname,fi)
            if sect.get(sctname) != None:
                sect[sctname].append(subsect)
            else:
                sect[sctname] = [subsect]
            continue
        key,val = line.split(maxsplit=1)
        if sect.get(key) != None:
            sect[key].append(val)
        else:
            sect[key] = [val]

def write_cp2k(sect_dict):
    lines = []
    def write_rec(sect_dict,l):
        for f in sect_dict:
            for i in sect_dict[f]:
                if type(i) == dict:
                    lines.append(' '*l+'&'+f)
                    l+=2
                    write_rec(i,l)
                    l-=2
                    lines.append(' '*l+'&END '+f)
                else:
                    if f == '':
                        lines.append(''*l+i.strip())
                    else:
                        lines.append(' '*l+f+' '+i)
        return '\n'.join(lines)+'\n'
    return(write_rec(sect_dict,1))

def write_input_cp2k(parm,fname='cp2k.inp'):
    with open(fname,'w') as fo:
        fo.write(write_cp2k(parm))

# import json

# sections = read_cp2k('cp2k.inp')
# print(sections)

# print(sections['FORCE_EVAL'][0]['MM'][0]['FORCEFIELD'][0]['FORCE_SCALE'][0])
# print(json.dumps(sections,sort_keys=True, indent=4))

# exit()

# with open('cp2k2.inp','w') as fo:
#     txt = write_cp2k(sections)
#     fo.write(txt)






def reverse_readline(filename, buf_size=8192):
    """A generator that returns the lines of a file in reverse order"""
    with open(filename) as fh:
        segment = None
        offset = 0
        fh.seek(0, os.SEEK_END)
        file_size = remaining_size = fh.tell()
        while remaining_size > 0:
            offset = min(file_size, offset + buf_size)
            fh.seek(file_size - offset)
            buffer = fh.read(min(remaining_size, buf_size))
            remaining_size -= buf_size
            lines = buffer.split('\n')
            # The first line of the buffer is probably not a complete line so
            # we'll save it and append it to the last line of the next buffer
            # we read
            if segment is not None:
                # If the previous chunk starts right from the beginning of line
                # do not concat the segment to the last line of new chunk.
                # Instead, yield the segment first
                if buffer[-1] != '\n':
                    lines[-1] += segment
                else:
                    yield segment
            segment = lines[0]
            for index in range(len(lines) - 1, 0, -1):
                if lines[index]:
                    yield lines[index]
        # Don't yield None if the file was empty
        if segment is not None:
            yield segment

def update_parm(new,base=DEFAULT_PARM):
    sections = new.keys() & base.keys()
    assert len(sections.difference(new.keys()))==0, f"unknown section(s) in: {new}"
    for section in sections:
        base[section].update(new[section])
    return base

def write_input(parm,fname='xtb.inp'):
    with open(fname,'w') as fo:
        for section in parm:
            fo.write(section+'\n')
            subsect = parm[section]
            for key in subsect:
                sign = SPECIAL_SIGNS.get(key,'=')
                fo.write(f'  {key}{sign}{subsect[key]}\n')
        fo.write('$end')

def parse_input(fname='xtb.inp'):
    parm = {}
    keys = {}
    section = None
    with open(fname,'r') as fi:
        lines = fi.readlines()
    delim = set(SPECIAL_SIGNS.values())
    delim.add('=')
    pat = f'[{"".join(delim)}]'
    # pat = re.compile(f'\s+([a-z]+)\s*[{",".join(delim)}]\s*(.*)')
    for line in lines:
        if '$' in line:
            if section!=None:
                parm[section] = keys.copy()
            section = line.strip()
            keys.clear()
        else:
            kw,val = re.split(pat,line.strip())
            keys[kw] = val
    return parm

def get_frame_xyz(fname,index=-1,frames=None,as_np=False,nmax=1000):
    xyz = []
    na = 0
    if index == -1:
        lines = iter(reverse_readline(fname))
        while True:
            line = next(lines).strip()
            # line = line.strip()
            if 'xtb' in line or 'time' in line or 'E =' in line:
                comment = line
                break
            if na > nmax:
                raise RuntimeError(f'incorrect .xyz file: {fname}, check comment line \
                    (must contain "xtb") and N atoms (if N > 1k, add "nmax = N" to call)')
            xyz.append(line)
            na+=1
        na0 = next(lines).strip()
        assert int(na0)==na, f'incorrect N atoms: {na} was read, expected {na0}'
        if as_np:
            xyz.reverse()
            xyz_np = np.array([l.split()[1:] for l in xyz],dtype=float)
            return xyz_np
        else:
            xyz.append(comment)
            xyz.append(na0)
            xyz.reverse()
            return '\n'.join(xyz)
    else:
        assert index>=0
        with open(fname,'r') as fi:
            nframes = 0
            eof = False
            line1 = fi.readline().strip()
            nat = int(line1)
            buf = collections.deque([line1],nat+3)
            while True:
                line = fi.readline().strip()
                if 'xtb' in line or 'time' in line:
                    nframes+=1
                if nframes > index+1:
                    break
                if line == '':
                    eof = True
                    break
                buf.append(line.strip())
            assert nframes > 0, f'incorrect .xyz file: {fname}, check comment line (must contain "xtb")'
            if eof:
                buf.popleft()
            else:
                buf.pop()
            if as_np:
                buf.popleft() # - natoms
                buf.popleft() # - comment
                xyz_np = np.array([l.split()[1:] for l in buf],dtype=float)
                return xyz_np
            else:
                return '\n'.join(buf)

def get_frame_xyz2(fname,index=-1,as_np=False): # TODO replace get_frame_xyz
    nread = 0
    xyz = []
    assert type(index)==int or type(index)==list, 'provide integer or list'
    if type(index) == list:
        return
    if index < 0:
        print('not implemented')
        return
    with open(fname,'r') as fi:
        while nread < index:
            lines = []
            try:
                nat = int(fi.readline())
            except:
                # if nread==0:
                #     print('EOF reached or incorrect number of atoms!')
                # else:
                #     print(f'Trajectory with {n} structures was read!')
                break
            comm = fi.readline()
            lines.append(f'{nat}\n')
            lines.append(comm)
            for i in range(nat):
                lines.append(fi.readline())
            struc = ''.join(lines)
            # structures.append(struc)
            n+=1

    if index == -1:
        lines = iter(reverse_readline(fname))
        while True:
            line = next(lines).strip()
            # line = line.strip()
            if 'xtb' in line or 'time' in line:
                comment = line
                break
            if na > nmax:
                raise RuntimeError(f'incorrect .xyz file: {fname}, check comment line \
                    (must contain "xtb") and N atoms (if N > 1k, add "nmax = N" to call)')
            xyz.append(line)
            na+=1
        na0 = next(lines).strip()
        assert int(na0)==na, f'incorrect N atoms: {na} was read, expected {na0}'
        if as_np:
            xyz.reverse()
            xyz_np = np.array([l.split()[1:] for l in xyz],dtype=float)
            return xyz_np
        else:
            xyz.append(comment)
            xyz.append(na0)
            xyz.reverse()
            return '\n'.join(xyz)
    else:
        assert index>=0
        with open(fname,'r') as fi:
            nframes = 0
            eof = False
            line1 = fi.readline().strip()
            nat = int(line1)
            buf = collections.deque([line1],nat+3)
            while True:
                line = fi.readline().strip()
                if 'xtb' in line or 'time' in line:
                    nframes+=1
                if nframes > index+1:
                    break
                if line == '':
                    eof = True
                    break
                buf.append(line.strip())
            assert nframes > 0, f'incorrect .xyz file: {fname}, check comment line (must contain "xtb")'
            if eof:
                buf.popleft()
            else:
                buf.pop()
            if as_np:
                buf.popleft() # - natoms
                buf.popleft() # - comment
                xyz_np = np.array([l.split()[1:] for l in buf],dtype=float)
                return xyz_np
            else:
                return '\n'.join(buf)

# xyz = get_frame_xyz('example/xtb.trj',-1)
# with open('test-m.xyz','w') as fo:
#     fo.write(xyz)
#     fo.write('\n')
#     fo.write(xyz)

# write_input(DEFAULT_PARM,'xtb.inp')
# d = parse_input()
# write_input(d,'xtb2.inp')

def process_traj_xyz(fname, traj_pdb=None, dump_pdb=False, dump_sep_cif=False,cell=None,ref_cif=None):
    # if traj_pdb==None:
    #     traj_pdb = fname[:-4]+'.pdb'
    # open(traj_pdb,'w').close()
    with open(fname,'r') as fi:
        nframes = 0
        if dump_sep_cif:
            try:
                os.mkdir('cif')
            except:
                shutil.rmtree('cif')
                os.mkdir('cif')
        while True:
            nread = 0
            line1 = fi.readline().strip()
            if line1 == '':
                break
            nat = int(line1)
            if nframes==0:
                buf = collections.deque([],nat)
            line2 = fi.readline().strip()
            pat = 'E =\s*(-*[0-9]*\.[0-9]+)'
            if re.search(pat,line2):
                ener = re.findall(pat,line2)[0]
            else:
                ener = 'NaN'
            while nread<nat:
                line = fi.readline().strip()
                buf.append(line.strip())
                nread+=1
            nframes+=1
            xyz = np.array([line.split()[1:] for line in buf],dtype=float)
            sym = np.array([line.split()[1] for line in buf],dtype=float)
            if dump_sep_cif:
                assert cell!=None, 'could not create .cif without cell parameters!'
                assert ref_cif!=None, 'could not create .cif without reference!'

                M = get_cartesian_to_fractional_matrix(
                    cell['a'],
                    cell['b'],
                    cell['c'],
                    cell['alpha'],
                    cell['beta'],
                    cell['gamma']
                    )
                fxyz = np.dot(xyz,M)
                os.chdir('cif')
                write_cif(f'{os.path.basename(fname)[:-4]}-{nframes}.cif',fxyz,ref_cif,nat)
                os.chdir('../')
            if dump_pdb:
                pass

def write_cif(fname, xyz, ref_cif, natom):
    xyz_part = False
    n = 0
    with open(ref_cif,'r') as fi:
        ref = fi.readlines()
    with open(fname,'w') as fo:
        for line in ref:
            if '_atom_site_fract_z' in line:
                xyz_part = True
                fo.write(line)
                continue
            if xyz_part and n<natom:
                lab,sym,_,_,_ = line.split()
                x,y,z = xyz[n]
                line = f'{lab} {sym} {x:.6f} {y:.6f} {z:.6f}\n'
                n+=1
            fo.write(line)


# os.chdir('/home/artem/PAM/20.10-1ns/')


# process_traj_xyz(
#     'structures.xyz',
#     traj_pdb='structures.pdb',
#     dump_pdb=True,
#     dump_sep_cif=True,
#     cell={
#         'a':11.8050,
#         'b':17.1640,
#         'c':7.3930,
#         'alpha':90.00,
#         'beta':90.00,
#         'gamma':90.00
#     },
#     ref_cif='/home/artem/HREMD/src/paracetamol-orthorhombic-P1.cif'
# )

# print(0)



def xyz_to_pdb(traj_xyz, cell, traj_pdb=None, as_np=False):
    if traj_pdb==None:
        traj_pdb = traj_xyz[:-4]+'.pdb'
    open(traj_pdb,'w').close()
    with open(traj_xyz,'r') as fi:
        nframes = 0
        while True:
            nread = 0
            line1 = fi.readline().strip()
            if line1 == '':
                break
            nat = int(line1)
            if nframes==0:
                buf = collections.deque([],nat)
            line2 = fi.readline().strip()
            pat = 'E =\s*(-*[0-9]*\.[0-9]+)'
            if re.search(pat,line2):
                ener = re.findall(pat,line2)[0]
            else:
                ener = 'NaN'
            while nread<nat:
                line = fi.readline().strip()
                buf.append(line.strip())
                nread+=1
            nframes+=1
            with open(traj_pdb,'a') as fo:
                fo.write(f'REMARK    Step {nframes}, E = {ener}\n')
                fo.write(f'CRYST1   {cell["a"]:.3f}   {cell["b"]:.3f}    {cell["c"]:.3f}\
                           {cell["alpha"]:.2f}  {cell["beta"]:.2f}  {cell["gamma"]:.2f}\n')
                nwrite = 0
                for line in buf:
                    nwrite+=1
                    spl = line.split()
                    sym = spl.pop(0)
                    x,y,z=[float(i) for i in spl]
                    pdb_line = 'ATOM'+str(nwrite).rjust(7)+' '+sym.ljust(4)+f'{x:.3f}'.rjust(22)+f'{y:.3f}'.rjust(8)+f'{z:.3f}'.rjust(8)+'  0.00  0.00           '+sym+'\n'
                    fo.write(pdb_line)
                fo.write('END\n')

# os.chdir('/home/artem/paracetamol/pdb_vs_xyz/')
# xyz_to_pdb('traj.xyz',traj_pdb='test.pdb')

CELL={
            'a':11.8050,
            'b':17.1640,
            'c':7.3930,
            'alpha':90.00,
            'beta':90.00,
            'gamma':90.00,
            'SpGr': 'Pcab'
        }


# os.chdir('/home/artem/PAM/24.10-6ns/')
xyz_to_pdb('/home/artem/PAM/24.10-6ns/1.000e+00/min.xyz',CELL)
# print(0)

if TEST:
    import time
    #xyz (eof)
    start = time.time()
    xyz = get_frame_xyz('example/xtb.trj',-1)
    end = time.time()
    print(f'from EOF: {end-start} seconds')
                
    #xyz (bof)
    start = time.time()
    xyz = get_frame_xyz('example/xtb.trj',999)
    end = time.time()
    print(f'from BOF: {end-start} seconds')

    #xyz numpy
    xyz_np = get_frame_xyz('example/xtb.trj',as_np=True)
    print(xyz_np)

    # write/parse
    write_input(DEFAULT_PARM,'xtb.inp')
    d = parse_input()
    write_input(d,'xtb2.inp')