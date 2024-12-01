import os, shutil, re, subprocess, time
import base_utils

LAMMPS_PATH = '/home/artem/LAMMPS/lammps-static/bin/'

PLUMED_INPUT = '''m: DIFF ...
	
	PATTERN_FILE=%(xray_data)s
	NAMES_FILE=%(atom_names)s
	FORCE_TYPE=%(force_type)s
	FORCE_COEFF=%(force_const)s

	CELL_A=%(cell_a)s
	CELL_B=%(cell_b)s
	CELL_C=%(cell_c)s
		
	CELL_ALPHA=%(cell_alpha)s 
	CELL_BETA=%(cell_beta)s
	CELL_GAMMA=%(cell_gamma)s
        GAP=%(gap)s

...


# COLVAR - file with xray diff
PRINT ARG=m STRIDE=1 FILE=COLVAR
'''


VARS = {
    'N' : 2500,
    'alpha' : 1.0,
    'S' : 999999,
    'T' : 273.15,
    'd' : '',
    'r' : '',
    't' : '',
    'c' : '',
    'ndump' : 0,
}

KEYS = {
    'units' : 'real',
    'atom_style' : 'full',
    'dimension' : '3',
    'boundary' : 'p p p',
    'special_bonds' : 'lj 0.0 0.0 0.5 coul 0.0 0.0 0.8333333333',
    'pair_style' : 'lj/cut/coul/long/soft 1.0 0.5 10.0 8.0 8.0',
    'bond_style' : 'harmonic',
    'angle_style' : 'harmonic',
    'dihedral_style' : 'opls',
    'improper_style' : 'cvff',
    'kspace_style' : 'pppm 1.0e-4',
    'kspace_modify' : 'gewald 0.001 compute no',
    'neighbor' : '2.0 bin',
    'neigh_modify' : 'delay 10 every 1 check yes once no',
    'timestep' : '2.0',
    'velocity' : 'all create $T $S dist gaussian mom no rot no',
    'thermo_style' : 'custom time temp pe ke press vol density',
    'thermo' : '10',
    'fix 1' : 'all nvt temp $T $T 20.0',
    'fix 2' : 'all adapt 0 pair lj/cut/coul/long/soft lambda * * v_alpha reset yes scale yes',
    'dump' : 'DUMPFILE all xyz 1 $t',
    'write_dump' : ' all xyz $t modify element C N C C H C H C O H C H C H H C H H H O',
    'dump_modify' : 'DUMPFILE element C N C C H C H C O H C H C H H C H H H O',
    'min_style' : 'cg',
    'min_modify' : 'dmax 0.2',
    'minimize' : '1.0e-4 1.0e-6 500 1000'
}

COMP = {
    'compute peratom' : 'all pe/atom'
}


THERMO_PATT = {
    'time':'Time',
    'temp':'Temp',
    'pe':'PotEng',
    'ke':'KinEng',
    'press':'Press',
    'vol':'Volume',
    'density':'Density'
}

    # TODO implement energy minimization
def runMD(args):
    WD, IN_FNAME, DATA_FNAME, ERR_FNAME, TRJ_FNAME, K, PATT, ncores,  *_ = args
    os.chdir(WD)
    assert os.path.exists(IN_FNAME), 'Input file not found'
    assert os.path.exists(DATA_FNAME), 'Data file not found'
    # CMD = [os.path.join(LAMMPS_PATH,'lmp'), '-i', IN_FNAME]
    # CMD = ['apptainer', 'exec', '/home/users/igorstan/lammps_plumed.sif', '/opt/lammps/build/lmp', '-i', IN_FNAME]
    CMD = ['lmp', '-i', IN_FNAME]
    if ncores > 1:
        pass
        # CMD = ['mpirun','-n',f'{ncores}','lmp', '-i', IN_FNAME]
    thl = False
    nstep = 0
    dE, dV = [], []
    t0 = time.time()
    proc = subprocess.Popen(CMD, stdout=subprocess.PIPE, universal_newlines=True)
    # print(f'WD: {WD}; NP={ncores}')
    for stdout_line in iter(proc.stdout.readline, ""):
        # print(stdout_line)
        if re.search(PATT,stdout_line) != None:
            thl = True
            continue
        if re.search('Performance',stdout_line) !=None and thl: # TODO: remove?
            break
        if thl:
            try:
                # TODO: extend for unknown number and order of columns
                _, _, e, _, _, _, _ = [int(stdout_line.split()[0])] + [float(i) for i in stdout_line.split()[1:]]
                bias = float('NaN')
            except:
                continue
            nstep+=1
            dE.append(e)
            dV.append(bias)
    exc = proc.wait()
    proc.stdout.close()
    t = time.time()
    if exc:
        with open(ERR_FNAME,'a') as err:
            err.write(f'\ncalculation crashed at {nstep} step\n\
                        LAMMPS command line: {" ".join(CMD)}\n')
            with open(IN_FNAME,'r') as fi:
                err.write(f'input file: \n\n{fi.read()}\n')
        minEStruc, lastStruc, results = {}, {}, {}
        minEStruc['energy'], minEStruc['bias'], minEStruc['xyz'] = None, None, None
        lastStruc['energy'], lastStruc['bias'], lastStruc['xyz'] = None, None, None
        results['time'], results['exitCode'], results['eShift'], \
            results['biasShift'], results['lastStruc'], results['minEStruc'] = \
                None, exc, None, None, lastStruc, minEStruc
        return results
    else:
        assert os.path.exists(TRJ_FNAME), 'Trajectory file not found'
        # nMin = min(range(len(dE)), key=dE.__getitem__)
        nMin = min(range(len(dE)//2,len(dE)), key=dE.__getitem__)
        # print(dE)
        minEStruc, lastStruc, results = {}, {}, {}
        minEStruc['energy'], minEStruc['bias'], minEStruc['xyz'] = \
            dE[nMin], dV[nMin], base_utils.readXYZ(TRJ_FNAME,nMin*K)
        lastStruc['energy'], lastStruc['bias'], lastStruc['xyz'] = \
            dE[-1], dV[-1], base_utils.readXYZ(TRJ_FNAME, -1)
        results['time'], results['exitCode'], results['eShift'], \
            results['biasShift'], results['lastStruc'], results['minEStruc'] = \
                t-t0, exc, dE, dV, lastStruc, minEStruc
        return results

def run_single_point(args) -> float:
    WD, IN_FNAME, DATA_FNAME, ERR_FNAME, TRJ_FNAME, K, PATT, ncores, *_ = args
    os.chdir(WD)
    assert os.path.exists(IN_FNAME), 'Input file not found'
    assert os.path.exists(DATA_FNAME), 'Data file not found'
    # CMD = [os.path.join(LAMMPS_PATH,'lmp'), '-i', IN_FNAME]
    CMD = ['lmp', '-i', IN_FNAME]
    # CMD = ['/opt/lammps/build/lmp', '-i', IN_FNAME]
    # CMD = ['apptainer', 'exec', '/home/users/igorstan/lammps_plumed.sif', '/opt/lammps/build/lmp', '-i', IN_FNAME]
    if ncores > 1:
        pass
        # CMD = ['mpirun','-n',f'{ncores}','lmp', '-i', IN_FNAME]
        # print('MPI parallelization not implemented (k4)')
        # CMD = ['mpirun','-n',f'{ncores}','apptainer exec /home/users/igorstan/lammps_plumed.sif /opt/lammps/build/lmp', '-i', IN_FNAME]
    t0 = time.time()
    e = float('NaN')
    proc = subprocess.Popen(CMD, stdout=subprocess.PIPE, universal_newlines=True)
    iterstdout = iter(proc.stdout.readline, "")
    thl = False
    # print(f'\n\n\nWD: {WD}\n\n')
    for stdout_line in iterstdout:
        # print(stdout_line)
        if re.search(PATT,stdout_line) != None:
            thl = True
            continue
        if thl:
            try:
                # TODO: extend for unknown number and order of columns
                _, _, e, _, _, _, _ = [int(stdout_line.split()[0])] + [float(i) for i in stdout_line.split()[1:]]
                break
            except:
                print('Error in pattern of thermo_style header')
                break
    exc = proc.wait()
    proc.stdout.close()
    t = time.time()
    if exc:
        print('WARNING: single point calculation failed for some reason')
    return e


def runOpt(args):
    WD, IN_FNAME, DATA_FNAME, ERR_FNAME, TRJ_FNAME, K, PATT, ncores, *_ = args
    os.chdir(WD)
    assert os.path.exists(IN_FNAME), 'Input file not found'
    assert os.path.exists(DATA_FNAME), 'Data file not found'
    # CMD = [os.path.join(LAMMPS_PATH,'lmp'), '-i', IN_FNAME]
    CMD = ['lmp', '-i', IN_FNAME]
    # CMD = ['apptainer', 'exec', '/home/users/igorstan/lammps_plumed.sif', '/opt/lammps/build/lmp', '-i', IN_FNAME]
    if ncores > 1:
        pass
        # CMD = ['mpirun','-n',f'{ncores}','lmp', '-i', IN_FNAME]
    t0 = time.time()
    proc = subprocess.Popen(CMD, stdout=subprocess.PIPE, universal_newlines=True)
    iterstdout = iter(proc.stdout.readline, "")
    for stdout_line in iterstdout:
        if re.search('Energy initial, next-to-last, final',stdout_line) !=None:
            final = next(iterstdout)
            e0, _,  e1 = [float(i) for i in final.split()]
    exc = proc.wait()
    proc.stdout.close()
    t = time.time()
    if TRJ_FNAME=='':
        return
    if exc:
        e0 = e1 = float('NaN')
        xyz = None
        print('WARNING: optimization failed for some reason')
    else:
        xyz = base_utils.readXYZ(TRJ_FNAME,0)
    return (e0,e1,xyz)
    
    




class Simulation:
    IN_FNAME = 'lammps.inp'
    RST_FNAME = 'lammps.restart'
    TRJ_FNAME = 'traj.xyz'
    ERR_FNAME = 'lammps.error'
    DAT_FNAME = 'lammps.data'
    PLUMED_FNAME = 'plumed.inp'

    OPTS = {
        'minBeforeMD' : True,
        'minAfterMD' : False,
        'doLongES' : True,
        'doNPT' : False,
        'unwrapXYZ' : False,
        'optimize' : False,
        'compute PE' : False,
        'blank' : False,
        'single_point' : False,
        'mtd_xrd' : False,
        'patch_xyz' : False # to solve LAMMPS problem with labelmap when parsing xyz file using `read_dump` command
    }
    def __init__(self,alpha,WD,datfile,xyz = None, ndump = None, parm = None, vars = None, options = None,**kwargs):

        self.restart = False
        self._nrst = 0
        self.patt = None
        self._n = None # total number of MD steps
        self._nth = None # number of properties calculations in trajectory = _n/_nth
        self._ntr = None # number of frames in trajectory = _n/_ntr
        self.k = None # self._nth//self._ntr
        # TODO numpy arrays with values from each MD step
        self.ener = None
        self.bias = None
        self.rmsd = None
        self.dxrd = None
        self.results = None # dict with MD information (TODO: replace with separate class)
        
        assert os.path.isfile(datfile), f'Data file not found: {datfile}'
        
        # PLUMED data
        self.plumed = kwargs.get('plumed',{})
        # cell information
        self.cell_a = kwargs.get('cell_a',None)
        self.cell_b = kwargs.get('cell_b',None)
        self.cell_c = kwargs.get('cell_c',None)
        self.cell_alpha = kwargs.get('cell_alpha',None)
        self.cell_beta = kwargs.get('cell_beta',None)
        self.cell_gamma = kwargs.get('cell_gamma',None)
        # parameters (below)
        
        
        if None in [self.cell_a, self.cell_b, self.cell_c, self.cell_alpha, self.cell_beta, self.cell_gamma]:
            print('read unit cell from datafile')
            self.cell_a, self.cell_b, self.cell_c, self.cell_alpha, self.cell_beta, self.cell_gamma = base_utils.get_cell_lammps(datfile)
        else:
            print('using existing unit cell parameters')
        
        self.plumed.update({'cell_a':self.cell_a,'cell_b':self.cell_b,'cell_c':self.cell_c,'cell_alpha':self.cell_alpha,'cell_beta':self.cell_beta,'cell_gamma':self.cell_gamma})
        
        try:
            with open(datfile,'r') as f:
                f.read()
            print(f'Running simulation from scratch: {datfile}')
        except UnicodeDecodeError:
            print(f'Restarting simulation from file: {datfile}')
            self.restart = True

        if os.path.exists(os.path.join(WD,self.RST_FNAME)):
            self._nrst+=1
            print('Previous restart files detected, overwriting avoided')
            while True:
                rname = os.path.join(WD,self.RST_FNAME+f'-{self._nrst}')
                if os.path.exists(rname):
                    self._nrst+=1
                else:
                     break
        else:
            rname = os.path.join(WD,self.RST_FNAME)
        try:
            os.chdir(WD)
        except FileNotFoundError:
            print(f'Creating directory from scratch: {WD}')
            os.mkdir(WD)
        self.alpha = alpha
        self.WD = WD
        self.datfile = datfile
        self.vars = VARS.copy()
        self.keys = KEYS.copy()
        self.options = self.OPTS.copy()
        if parm != None:
            self._updateParm(parm)
            self._updateParm(kwargs.get('parm',{}))
        if vars != None:
            self._updateVar(vars)
        if options != None:
            self._updateOpt(options) #FIXME - only kwargs
            self._updateParm(kwargs.get('options',{}))
        self.vars['d'] = datfile
        self.vars['r'] = rname
        self.vars['t'] = self.TRJ_FNAME
        self.vars['alpha'] = alpha
        if xyz != None:
            self.vars['c'] = xyz
        if ndump != None:
            self.vars['ndump'] = ndump

    def _updateParm(self,parm):
        assert type(parm) is dict, 'dictionary object expected here'
        assert set(parm).issubset(set(self.keys)), f'Unknown keywords provided: {", ".join([f"{i}" for i in set(parm)-set(self.keys)])}'
        self.keys.update(parm)
        print('Default parameters were updated')

    def _updateVar(self,vars):
        assert type(vars) is dict, 'dictionary object expected here'
        assert set(vars).issubset(set(self.vars)), f'Unknown variables provided: {", ".join([f"{i}" for i in set(vars)-set(self.vars)])}'
        self.vars.update(vars)
        print('Variables were updated')

    def _updateOpt(self,options):
        assert type(options) is dict, 'dictionary object expected here'
        assert set(options).issubset(set(self.options)), f'Unknown options provided: {", ".join([f"{i}" for i in set(options)-set(self.options)])}'
        self.options.update(options)
        print('LAMMPS options were updated')


    def _checkParm(self):
        self._n = int(self.vars['N'])
        self._nth = int(self.keys['thermo'])
        self._ntr = int(self.keys['dump'].split()[3])
        thermst = self.keys['thermo_style'].split()
        assert 'custom' == thermst[0], 'custom thermo_style required. For example: "thermo_style     custom time temp pe ke press vol density"'
        assert set(thermst[1:]).issubset(set(THERMO_PATT)), f'unknown properties in thermo_style: {", ".join([f"{i}" for i in set(thermst[1:])-set(THERMO_PATT)])}'
        assert 'pe' in thermst, 'potential energy declaration in thermo_style required'
        self.patt = re.compile('\s+'.join([THERMO_PATT[p] for p in thermst[1:]]))
        self.k = self._nth//self._ntr


    # TODO implement energy minimization
    def _runLAMMPS(self,queue = None, structure = None):
        os.chdir(self.WD)
        assert os.path.exists(self.IN_FNAME), 'Input file not found'
        assert os.path.exists(self.vars['d']), 'Data file not found'
        # CMD = [os.path.join(LAMMPS_PATH,'lmp'), '-i', self.IN_FNAME]
        CMD = ['lmp', '-i', self.IN_FNAME]
        thl = False
        nstep = 0
        dE, dV = [], []
        t0 = time.time()
        proc = subprocess.Popen(CMD, stdout=subprocess.PIPE, universal_newlines=True)
        for stdout_line in iter(proc.stdout.readline, ""):
            if re.search(self.patt,stdout_line) != None:
                thl = True
                continue
            if re.search('Performance',stdout_line) !=None and thl: # TODO: remove?
                break
            if thl:
                try:
                    # TODO: extend for unknown number of columns and their order
                    _, _, e, _, _, _, _ = [int(stdout_line.split()[0])] + [float(i) for i in stdout_line.split()[1:]]
                    bias = float('NaN')
                except:
                    continue
                nstep+=1
                dE.append(e)
                dV.append(bias)
        exc = proc.wait()
        proc.stdout.close()
        t = time.time()
        if exc:
            with open(self.ERR_FNAME,'a') as err:
                err.write(f'\ncalculation crashed at {nstep} step\n\
                          LAMMPS command line: {" ".join(CMD)}\n')
                with open(self.IN_FNAME,'r') as fi:
                    err.write(f'input file: \n\n{fi.read()}\n')
            minEStruc, lastStruc, results = {}, {}, {}
            minEStruc['energy'], minEStruc['bias'], minEStruc['xyz'] = None, None, None
            lastStruc['energy'], lastStruc['bias'], lastStruc['xyz'] = None, None, None
            results['time'], results['exitCode'], results['eShift'], \
                results['biasShift'], results['lastStruc'], results['minEStruc'] = \
                    None, exc, None, None, lastStruc, minEStruc
            return results
        else:
            assert os.path.exists(self.TRJ_FNAME), 'Trajectory file not found'
            nMin = min(range(len(dE)), key=dE.__getitem__)
            minEStruc, lastStruc, results = {}, {}, {}
            minEStruc['energy'], minEStruc['bias'], minEStruc['xyz'] = \
                dE[nMin], dV[nMin], base_utils.readXYZ(self.TRJ_FNAME,nMin*self.k)
            lastStruc['energy'], lastStruc['bias'], lastStruc['xyz'] = \
                dE[-1], dV[-1], base_utils.readXYZ(self.TRJ_FNAME, -1)
            results['time'], results['exitCode'], results['eShift'], \
                results['biasShift'], results['lastStruc'], results['minEStruc'] = \
                    t-t0, exc, dE, dV, lastStruc, minEStruc
            return results

    def _writeInput(self):
        block1 = ['units', 'atom_style', 'timestep', 'dimension', 'boundary', 'special_bonds', 'pair_style', 'bond_style', 'angle_style', 'dihedral_style', 'improper_style']
        block2 = ['velocity','kspace_style', 'kspace_modify', 'neighbor', 'neigh_modify']
        block3 = ['dump', 'dump_modify', 'thermo_style', 'thermo', 'fix 1', 'fix 2']
        minimize = ['min_style', 'min_modify', 'minimize']
        unused_opt = ['velocity', 'neighbor', 'neigh_modify', 'thermo_style', 'thermo', 'fix 1', 'fix 2', 'timestep']
        unused_restart = ['velocity'] + block1
        os.chdir(self.WD)
        keys = self.keys.copy()
        # apply options
        if self.options['doLongES']:
            keys['kspace_modify'] = keys['kspace_modify'].replace('compute no','compute yes')
        else:
            keys['kspace_modify'] = keys['kspace_modify'].replace('compute yes','compute no')
        if self.options['doNPT']:
            keys['fix 1'] = 'all npt temp $T $T 100 iso 1 1 1000'
        if self.options['unwrapXYZ']:
            keys['dump'] = 'DUMPFILE all custom 1 $t element xu yu zu'
            print('WARNING: atomic coordinates will be printed in unwrapped format')
        if self.options['mtd_xrd']:
            if self._check_plumed():
                print(f'XRD metadynamics requested. plumed file = {os.path.join(self.WD, self.PLUMED_FNAME)}, xray data: {self.plumed["xray_data"]}, output: {os.path.join(self.WD, self.PLUMED_FNAME[:-4]+".out")}')
                keys['fix mtd'] = f'all plumed plumedfile {os.path.join(self.WD, self.PLUMED_FNAME)} outfile {os.path.join(self.WD, self.PLUMED_FNAME[:-4]+".out")}'
                # fix mtd all plumed plumedfile example_plumed_input.dat outfile p.log
                block3.append('fix mtd')
                self._write_plumed()
            else:
                print('XRD metadynamics requested but not all parameters were provided, skip')
        comp = ''
        if self.options['compute PE']:
            comp = 'compute peratom ' + COMP['compute peratom'] + '\n'
        read = '\n\nread_data $d\n\n'
        if self.restart:
            block2 = [i for i in block2 if i not in unused_restart]
            read = '\nread_restart $d\n\n'
        if self.vars['c'] != '':
            read+=f"read_dump {self.vars['c']} {self.vars['ndump']} x y z box no format xyz\n\n"
            print(f"coordinates taken from {self.vars['c']}, frame {self.vars['ndump']}")
            self.vars['c'] = ''
            self.vars['ndump'] = 0
        if self.options['optimize']:
            block1 = [i for i in block1 if i not in unused_opt]
            block2 = [i for i in block2 if i not in unused_opt]
            minimize+=['write_dump']
        act_vars = [v for v in self.vars if v not in ['c','ndump']]
        _vars = [f'variable {v} string {self.vars[v]}\n' if type(self.vars[v]) is str \
                 else f'variable {v} equal {self.vars[v]}\n' for v in act_vars]
        rows1 = [f'{v} {keys[v]}\n'for v in block1]
        rows2 = [f'{v} {keys[v]}\n'for v in block2]
        rows3 = [f'{v} {keys[v]}\n'for v in block3]
        rowsm = [f'{v} {keys[v]}\n'for v in minimize]
        # build order
        if self.options['optimize']:
            order = _vars + ['\n\nclear\n\n'] + rows1 + [read] + rows2 + rowsm
        else:
            order = _vars + ['\n\nclear\n\n'] + rows1 + [read] + rows2 + [comp] + rows3 + ['\nreset_timestep 0\nrun $N\n\n\nwrite_restart $r']
            if self.options['minBeforeMD']:
                order = _vars + ['\n\nclear\n\n'] + rows1 + [read] + rows2 + [comp] + rowsm + rows3 + ['\nreset_timestep 0\nrun $N\n\n\nwrite_restart $r']
            if self.options['minAfterMD']:
                print('WARNING: minimization after MD not implemented, skip it')
        if self.options['blank']:
            print('Warning: no "run" or "minimize" task for lammps')
            order = _vars + ['\n\nclear\n\n'] + rows1 + [read] + rows2 + ['\n\n\nwrite_data %s\n\n' % os.path.join(self.WD, self.DAT_FNAME)]
        if self.options['single_point']:
            order = _vars + ['\n\nclear\n\n'] + rows1 + [read] + rows2 + rows3 + ['\nrun 0\n\n\nwrite_restart $r']            
        # write to input file
        with open(self.IN_FNAME,'w') as fi:
            for row in order:
                fi.write(row)

    # TODO: parallel implementation without division into 3 methods
    def prepare(self,newrestart = None, parm = None, vars = None, options = None):
        if newrestart != None:
            self.updateRestartFile(newrestart)
        if parm != None:
            self._updateParm(parm)
        if options != None:
            self._updateOpt(options)
        if vars != None:
            assert 'd' not in vars, 'uncertainty in the number of datafiles'
            self._updateVar(vars)
        self._checkParm()
        self._writeInput()

    def update(self,results):
        self._nrst+=1
        self.vars['d'] = self.vars['r']
        self.vars['r'] = os.path.join(self.WD, self.RST_FNAME+f'-{self._nrst}')
        self.results = results
        self.restart = True
        # self.options = self.OPTS.copy() # options are reset to defaults

    def _check_plumed(self):
        _required = set(['cell_alpha', 'cell_b', 'atom_names', 'cell_a', 'cell_beta', 'cell_gamma', 'force_const', 'cell_c', 'xray_data', 'force_type'])
        _diff = _required - set(self.plumed.keys())
        if len(_diff) > 0:
            print(f'following metadynamic parameters are missing: {" ".join(_diff)}')
            return False
        elif not os.path.exists(self.plumed['atom_names']):
            print(f'file with atom names not found: {self.plumed["atom_names"]}')
            return False
        elif not os.path.exists(self.plumed['xray_data']):
            print(f'diffraction data not found: {self.plumed["xray_data"]}')
            return False
        else:
            return True
        
        

    def _write_plumed(self):
        with open(os.path.join(self.WD,self.PLUMED_FNAME),'w') as fo:
            fo.write(PLUMED_INPUT % self.plumed)

    def runMD(self,newrestart = None, parm = None, vars = None, options = None):
        os.chdir(self.WD)
        if newrestart !=None:
            self.updateRestartFile(newrestart)
        if parm != None:
            self._updateParm(parm)
        if options != None:
            self._updateOpt(options)
        if vars != None:
            assert 'd' not in vars, 'uncertainty in the number of datafiles'
            self._updateVar(vars)
        self._checkParm()
        self._writeInput()
        self.results = self._runLAMMPS()
        self._nrst+=1
        self.vars['d'] = self.vars['r']
        self.vars['r'] = os.path.join(self.WD, self.RST_FNAME+f'-{self._nrst}')
        self.restart = True
        # self.options = self.OPTS.copy() # options are reset to defaults
    
    def updateRestartFile(self, newrestart):
        assert os.path.isfile(newrestart), f'Restart file not found: {newrestart}'
        self.vars['d'] = newrestart
        self.restart = True
        # print('restart file updated')

    def updateXyz(self, newxyz):
        with open(os.path.join(self.WD,'tmp.xyz'),'w') as fo:
            for row in newxyz:
                fo.write(row)
        if self.options['patch_xyz']:
            base_utils.patch_xyz(os.path.join(self.WD,'tmp.xyz'))
        self.vars['c'] = os.path.join(self.WD,'tmp.xyz')
        self.vars['ndump'] = 0


    def getRestartFile(self):
        rname = os.path.join(self.WD, self.vars['d'])
        assert self.restart == True, f'Restart file was requested before actual MD calculation: {rname}'
        assert os.path.exists(rname), f'Restart file was deleted or renamed: {rname}'
        return rname

    def _checkTraj(self):
        pass

    def clearDir(self,savelist=[]):
        exceptions = [os.path.basename(f) for f in [self.vars['d'], self.vars['r'], self.vars['t'], self.IN_FNAME, self.ERR_FNAME] + savelist]
        for file in os.listdir(self.WD):
            if file in exceptions:
                continue
            if os.path.isdir(file):
                try:
                    shutil.rmtree(file)
                except OSError as e:
                    print("Error: %s - %s." % (e.filename, e.strerror))
            else:
                try:
                    os.remove(file)
                except OSError as e:
                    print("Error: %s - %s." % (e.filename, e.strerror))

    def __del__(self):
        print(f'Simulation obj. destructed, dir: {self.WD}')
        # self.clearDir()

TEST = False
if TEST:
    opts = {
        'minBeforeMD':True,
        'doLongES':True
        }
    s1 = Simulation(0.0,'/home/artem/LAMMPS_TEST/dbg/0.0','/home/artem/LAMMPS_TEST/lammps.data')
    s2 = Simulation(0.5,'/home/artem/LAMMPS_TEST/dbg/0.5','/home/artem/LAMMPS_TEST/lammps.data', vars = {'N':600})
    s1.clearDir()
    s2.clearDir()
    s1.runMD(vars = {'N':600})
    s2.runMD()
    # s1.clearDir([])
    # s2.clearDir()
    r1 = s1.getRestartFile()
    r2 = s2.getRestartFile()
    s1.runMD(r2)
    s2.runMD(r1)
    r1 = s1.getRestartFile()
    r2 = s2.getRestartFile()
    s1.runMD(r2)
    s2.runMD(r1)

    print(0)

# TODO: CP2K interface
# TODO: XTB interface