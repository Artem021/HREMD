# low-level organisation of .xyz data
from collections.abc import Sequence
import random
import numpy as np
from io import StringIO
import os

import base_utils

# {coord[0]:8.3f} TODO

# super(Num2, self).__init__(num) TODO

class File(object):
    '''
    Representation of abstract text file (data or xyz or any but not binary/zip)
    use `super` when inheriting

    Arguments
    ---------
    `filename` : str
        Name of input file
    `mode` : str, optional
        read/write mode
    '''
    def __init__(self, filename: str, mode='r') -> None:
        self._filename = filename
        self._file = open(filename, mode)
        self.mode = mode
        self._is_open = True

    def guess_ext(self):
        return self._filename.split('.')[-1]
    
    def close(self) -> None:
        if self._is_open:
            self._file.close()
    
    def __enter__(self):
        return self
    
    def __exit__(self, *exc_) -> None:
        self.close()
    
    def __del__(self) -> None:
        self.close()
    
    def _validate(self) -> None:
        pass
    
    def open(self) -> None:
        pass
    

class XYZFile(File):
    
    def __init__(self, filename: str, mode='r') ->  None:
        super(XYZFile, self).__init__(filename, mode)
        self._file
        self._natoms = None
        self._nframes = None
        
    def read(self):
        pass
    
    @property
    def natoms(self):
        return self._natoms
    
    @property
    def nframes(self):
        return self._nframes


class LammpstrjFile(File):
    pass


def read_xyz(filename : str):
    frames = []
    _nf = 0
    with open(filename,'r') as fi:
        nlines = sum(1 for i in fi)
        fi.seek(0)
        line = fi.readline()
        if 'ITEM' in line:
            offset = 1
            while True:
                if 'ITEM: ATOMS' in line:
                    break
                if 'ITEM: NUMBER OF ATOMS' in line:
                    na = int(fi.readline())
                    offset+=1
                line = fi.readline()
                offset+=1
        else:
            offset = 2
            na = int(line)        
        assert nlines%(na+offset)==0, 'corrupted or nonstandart xyz file'
        nframes = nlines//(na+offset)
        fi.seek(0)
        lines = iter(fi.readlines())
        if index !=None: # иначе возвращаем все фреймы
            indices = []
            if type(index) is list:
                indices+=index
            else:
                indices.append(index)
            for i in range(len(indices)):
                ind = indices[i]
                assert ind<nframes, f'wrong index requested: {ind}'
                if ind<0:
                    ind+=nframes
                    assert ind>=0, f'wrong index requested: {indices[i]}'
                    indices[i]=ind
            indices = sorted(indices)
        else:
            indices = [i for i in range(nframes)]
        while True:
            for i in range(offset):
                try:
                    _ = next(lines)
                except StopIteration:
                    break
            if nf==indices[0]:
                xyz = []
                xyz.append(f'{na}\n')
                xyz.append(f'Frame {nf}\n')
                for i in range(na):
                    line = next(lines)
                    row = line.split()
                    el, x, y, z = row[-4:]
                    xyz.append(' '.join([el, x, y, z]) + '\n')
                    # pe = row[-1]
                    # xyz.append(' '.join([pe]) + '\n')
                frames.append(xyz)
                indices.pop(0)
                if len(indices)==0:
                    break
            else:
                for i in range(na):
                    _ = next(lines)
            nf+=1
    if len(frames)==1:
        return frames[0]
    return frames


def read_lammps():
    pass











class Default(object):
    def __getattr__(self, atrb) -> None:
        return None
    
# ?performance?
class Atom(Default):
    def __init__(self,x,y,z,sym) -> None:
        self.x = x
        self.y = y
        self.z = z
        self.sym = sym
    
    def translate(self,dx,dy,dz) -> None:
        self.x+=dx
        self.y+=dy
        self.z+=dz

    # dynamic properties may be arbitrary
    # ener, charge, alpha, atype, is_hboning, neigh, molid, sysid, image_flags

# 16432
# a_fl = np.array([[random.random(),random.random(),random.random()] for i in range(16432)])
#??? a_cust = np.array([Atom(random.random(),random.random(),random.random()) for i in range(16432)],dtype=Default)
# np.vectorize()

class Molecule(Default):
    # cell with PBC or not
    topo = None
    atoms = None
    
    def __init__(self, obj:list|str) -> None:
        super().__init__() # ???
        if isinstance(obj, list):
            print('creating new molecule from list')
            # self._from_list()
        elif isinstance(obj, str):
            print('creating new molecule from .xyz')
            # self._from_xyz()
        else:
            raise ValueError('unsupported data type for Molecule initialization, must be list or str')
        
    @staticmethod
    def from_file(cont : str) -> None:
        if os.path.isfile(cont):
            base_utils.readXYZ(cont,0)
        else:
            print('reading from string not implemented')
            raise NotImplementedError
        # _atoms = [Atom(), ...]
        pass
    
    def from_list(self):
        pass
    
    @property
    def index(self):
        return random.random()
    
    @property
    def xyz(self) -> np.ndarray:
        return np.array([])
        
    
    def wrapped(self) -> None:
        pass
    
    def unwrapped(self) -> None:
        pass
    
    def get_xyz(self) -> None:
        pass
    
    # @property
    def rmsd(self) -> float:
        print(f'rmsd is {random.random()}')
        return random.random()




class System(Default,Sequence): # base --> Cluster, ...??
    
    _atoms = None
    molecules = []
    
    @staticmethod
    def from_file() -> None:
        pass # Molecule
    
    def from_molecule(mol, n) -> None:
        pass # Molecule

    def from_molecules(mols, ns) -> None:
        pass # Molecule
    
    # lambda function --> set of atoms with desired properties
    def select(self): # returns some molecules/atoms
        pass
    
    def pack(self):
        pass
    
    def set_property(self,array, name): # set arbitrary property for each atom or molecule
        pass
    

# 
#
# 
# 
# 
# 
# 
# 
# 

# # Particle ?

# # 1 molecule repr

# xyz, connect_graph, cell, pbc, image_flags (ixyz)

# + map with indexes and symbols

# dynamic properties:

# dE (+contributions), distance M (closest contacts), arbitrary property

# dF

# wrap()

# unwrap()

# getxyz(wrapped=yes|no)

# get_rmsd(Structure())

# -------------------------------------

# # N molecules:

# # access by index??

# # common properties: unit cell

# System([Structure,'''])

# add()

# delete()

# calc_

# # high-level methods

# https://github.com/cvxgrp/dccp/tree/master

# pack() --> changes unit cell, struc.xyz, struc.ixyz



# идеи спизженные из mdtraj:

# 0. писать осмысленный код 

# 1. PDBTrajectoryFile - классы <-- файл с траекторией

# 2. 


# -------------------- expected -------------------- #
# growth_unit = Molecule()
# initial_cluster = Molecule()
# sys1 = System([initial_cluster, growth_unit])
# samples = get_samples(sys1)
# best = min(samples, lambda x: x.energy)
# sys2 = System(best, growth_init)


# reference rmsd
# get_ref_clusters(struc,n) -> list(Molecule)