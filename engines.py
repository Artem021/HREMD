import os, shutil, re, subprocess, time
import base_utils

LAMMPS_PATH = '/home/artem/LAMMPS/lammps-static/bin/'

VARS = {
    'N' : 2500,
    'alpha' : 1.0,
    'S' : 999999,
    'T' : 273.15,
    'd' : '',
    'r' : '',
    't' : ''
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
    'dump_modify' : 'DUMPFILE element C N C C H C H C O H C H C H H C H H H O',
    'min_style' : 'cg',
    'min_modify' : 'dmax 0.2',
    'minimize' : '1.0e-4 1.0e-6 500 1000'
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
def runLAMMPS(args):
    WD, IN_FNAME, DATA_FNAME, ERR_FNAME, TRJ_FNAME, K, PATT, *_ = args
    os.chdir(WD)
    assert os.path.exists(IN_FNAME), 'Input file not found'
    assert os.path.exists(DATA_FNAME), 'Data file not found'
    # CMD = [os.path.join(LAMMPS_PATH,'lmp'), '-i', IN_FNAME]
    CMD = ['lmp', '-i', IN_FNAME]
    thl = False
    nstep = 0
    dE, dV = [], []
    t0 = time.time()
    proc = subprocess.Popen(CMD, stdout=subprocess.PIPE, universal_newlines=True)
    for stdout_line in iter(proc.stdout.readline, ""):
        if re.search(PATT,stdout_line) != None:
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
        nMin = min(range(len(dE)), key=dE.__getitem__)
        minEStruc, lastStruc, results = {}, {}, {}
        minEStruc['energy'], minEStruc['bias'], minEStruc['xyz'] = \
            dE[nMin], dV[nMin], base_utils.get_frame_xyz(TRJ_FNAME,nMin*K)
        lastStruc['energy'], lastStruc['bias'], lastStruc['xyz'] = \
            dE[-1], dV[-1], base_utils.get_frame_xyz(TRJ_FNAME)
        results['time'], results['exitCode'], results['eShift'], \
            results['biasShift'], results['lastStruc'], results['minEStruc'] = \
                t-t0, exc, dE, dV, lastStruc, minEStruc
        return results


class Simulation:
    IN_FNAME = 'lammps.inp'
    RST_FNAME = 'lammps.restart'
    TRJ_FNAME = 'traj.xyz'
    ERR_FNAME = 'lammps.error'

    OPTS = {
        'minBeforeMD' : True,
        'minAfterMD' : False,
        'doLongES' : True
    }
    def __init__(self,alpha,WD,datfile,parm = None, vars = None, options = None):

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
        if vars != None:
            self._updateVar(vars)
        if options != None:
            self._updateOpt(options)
        self.vars['d'] = datfile
        self.vars['r'] = rname
        self.vars['t'] = self.TRJ_FNAME
        self.vars['alpha'] = alpha

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
        print('MD settings were updated')


    def _checkParm(self):
        self._n = int(self.vars['N'])
        self._nth = int(self.keys['thermo'])
        self._ntr = int(self.keys['dump'].split()[-2])
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
                dE[nMin], dV[nMin], base_utils.get_frame_xyz(self.TRJ_FNAME,nMin*self.k)
            lastStruc['energy'], lastStruc['bias'], lastStruc['xyz'] = \
                dE[-1], dV[-1], base_utils.get_frame_xyz(self.TRJ_FNAME)
            results['time'], results['exitCode'], results['eShift'], \
                results['biasShift'], results['lastStruc'], results['minEStruc'] = \
                    t-t0, exc, dE, dV, lastStruc, minEStruc
            return results

    def _writeInput(self):
        os.chdir(self.WD)
        keys = KEYS.copy()
        if self.options['doLongES']:
            keys['kspace_modify'] = 'gewald 0.001 compute yes'
        else:
            keys['kspace_modify'] = 'gewald 0.001 compute no'
        with open(self.IN_FNAME,'w') as fi:
            for v in self.vars:
                if type(self.vars[v]) is str:
                    eq = 'string'
                else:
                    eq = 'equal'
                fi.write(f'variable {v} {eq} {self.vars[v]}\n')
            if self.restart:
                fi.write('\nread_restart $d\n\n')
            else:
                fi.write('\nclear\n\n')
                for v in [
                    'units',
                    'atom_style',
                    'timestep',
                    'dimension',
                    'boundary',
                    'special_bonds',
                    'pair_style',
                    'bond_style',
                    'angle_style',
                    'dihedral_style',
                    'improper_style'
                    ]:
                    fi.write(f'{v} {keys[v]}\n')
                fi.write('\nread_data $d\n\n')
                fi.write(f'velocity {keys["velocity"]}\n')

            for v in [
                'kspace_style',
                'kspace_modify',
                'neighbor',
                'neigh_modify',
                'thermo_style',
                'thermo',
                'fix 1',
                'fix 2'
                ]:
                fi.write(f'{v} {keys[v]}\n')
            if self.options['minBeforeMD']:
                fi.write(f'min_style {keys["min_style"]}\n')
                fi.write(f'min_modify {keys["min_modify"]}\n')
                fi.write(f'minimize {keys["minimize"]}\n')
            for v in [
                'dump',
                'dump_modify'
                ]:
                fi.write(f'{v} {keys[v]}\n')

            fi.write('\nrun $N\n\n')
            if self.options['minAfterMD']:
                print('WARNING: minimization after MD not implemented, skip it')
            fi.write('\nwrite_restart $r')

    # TODO: parallel implementation except for dividing on 3 functions
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
    s2 = Simulation(0.5,'/home/artem/LAMMPS_TEST/dbg/0.5','/home/artem/LAMMPS_TEST/lammps.data', vars = {'N':500})
    s1.clearDir([])
    s2.clearDir(['ekhgeg897wg843'])
    s1.IN_FNAME = 'lammps.inp.0.0'
    s2.IN_FNAME = 'lammps.inp.0.5'
    s1.runMD(vars = {'N':600})
    s1.runMD()
    s1.clearDir([])
    s2.clearDir(['ekhgeg897wg843'])
    # s1.updateRestartFile('lammps.restart')
    # s1.runMD(options=opts,
    #     vars={
    #         'alpha':1.0, 
    #         'N':1000
    #         }
            # )
    s2.runMD()
    d = os.getcwd()
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