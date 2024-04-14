import re,os,subprocess,shutil,time,copy
import base_utils

XTB_PATH = '/home/artem/xtb-v2/build'
BLAS_PATH = '/home/artem/OpenBLAS/lib/'
LAMMPS_PATH = '/home/artem/LAMMPS/lammps-static/bin/'


class MD_calc:

    def __init__(self,DIR,XYZ,ARGS,CMD_LINE,INP,OUT,ERR,TRJ,RST,restart,alpha,Nxyz,XYZ_SHARED,engine) -> None:
        self.results = {
            'ener_shift' : [],
            'vb_shift' : [],
            'MD_steps' : 0,
            'error_code' : 0,
        }
        self.ener_shift = []
        self.vb_shift = []
        self.MD_steps = 0
        self.error_code = 0
        self.local_min_struc = {
            'Energy':None,
            'Vbias':None,
            'xyz':None
            }
        self.last_struc = {
            'Energy':None,
            'Vbias':None,
            'xyz':None
            }

        if not(os.path.exists(ERR)):
            open(ERR,'w').close()
        if not(os.path.exists(OUT)):
            open(OUT,'w').close()
        if engine == 'xtb':
            parm = base_utils.update_parm(ARGS)
            if Nxyz != 0:
                parm['$metadyn']['save'] = Nxyz
                if '--xyzset' not in CMD_LINE:
                    assert os.path.exists(XYZ_SHARED)
                    CMD_LINE.append('--xyzset')
                    CMD_LINE.append(XYZ_SHARED)
            parm.pop('$wall','')
            base_utils.write_input_xtb(parm,INP)
        elif engine == 'cp2k':
            try:
                _ = os.environ['OPENMPI_LIBS']
            except:
                raise EnvironmentError('type in console: source /home/artem/cp2k/tools/toolchain/install/setup')
            base_utils.write_input_cp2k(ARGS,fname=INP)
        # works only with template from /src/lammps.inp and /src/lammps.restart
        elif engine == 'lammps':
            if restart:
                assert RST != None
                base_utils.update_inp(os.path.join(DIR,RST),{'alpha':alpha})
            else:
                base_utils.update_inp(os.path.join(DIR,INP),{'alpha':alpha})
        self.DIR = DIR
        self.INP = INP
        self.OUT = OUT
        self.ERR = ERR
        self.XYZ = XYZ
        self.TRJ = TRJ
        self.RST = RST
        self.restart = restart
        self.ARGS = ARGS # dict with sections of input file
        self.CMD_LINE = CMD_LINE # [path/to/xtb, '--arg1', '--arg2', ...]
        self.alpha = alpha
        self.engine = engine

    def start_job(self,queue,index):
        if self.engine =='cp2k':
            return self.start_CP2K(queue,index)
        elif self.engine =='lammps':
            return self.start_LAMMPS(queue,index)
        elif self.engine =='xtb':
            return self.start_XTB(queue,index)
        
    def start_CP2K(self,queue,index):
        DIR = os.path.dirname(self.INP)
        os.chdir(DIR)
        pat = re.compile('ENERGY[=|\s]*(\S+)\s*VBIAS[=|\s]*(\S*)\s*MD_STEP[=|\s]*(\S*)')
        err = open(self.ERR,'a')
        t0 = time.time()
        popen = subprocess.Popen(self.CMD_LINE, stdout=subprocess.PIPE, universal_newlines=True,stderr=err)
        for stdout_line in iter(popen.stdout.readline, ""):
            res = re.search(pat,stdout_line)
            if res != None:
                e,vb,n = re.findall(pat,stdout_line)[0]
                try:
                    e = float(e)
                except:
                    e = float('NaN')
                try:
                    vb = float(vb)
                except:
                    vb = float('NaN')
                try:
                    n = int(n)
                except:
                    n = float('NaN')
                self.MD_steps+=1
                self.ener_shift.append(e)
                self.vb_shift.append(vb)
        self.return_code = popen.wait()
        popen.stdout.close()
        err.close()
        t = time.time()
        if self.return_code:
            with open(self.INP,'r') as fi:
                err_inp = fi.read()
            with open(self.XYZ,'r') as fi:
                err_xyz = fi.read()
            err = open(self.ERR,'a')
            err.write(f'\ncalculation is down at {self.MD_steps} step\n')
            err.write(f'command line args: {self.CMD_LINE}\n')
            err.write(f'input file: \n{err_inp}\n')
            err.write(f'xyz file: \n{err_xyz}\n')
            err.close()
        else:
            assert len(self.vb_shift) == len(self.ener_shift) == self.MD_steps
            # processing local minimum structure
            min_e = self.ener_shift[0]
            min_i = 0
            for i,e in enumerate(self.ener_shift):
                if e < min_e:
                    min_e = e
                    min_i = i
            min_vb = self.vb_shift[min_i]
            min_xyz = base_utils.get_frame_xyz(self.TRJ,min_i)
            self.local_min_struc['Energy'] = min_e
            self.local_min_struc['Vbias'] = min_vb
            self.local_min_struc['xyz'] = min_xyz
            self.local_min_struc['index'] = min_i
            self.local_min_struc['Nsteps'] = self.MD_steps
            # processing last structure
            last_e = self.ener_shift[-1]
            last_vb = self.vb_shift[-1]
            last_i = self.MD_steps-1
            last_xyz = base_utils.get_frame_xyz(self.TRJ)
            self.last_struc['Energy'] = last_e
            self.last_struc['Vbias'] = last_vb
            self.last_struc['xyz'] = last_xyz
            self.last_struc['index'] = last_i
        os.chdir('..')
        dt = t-t0
        queue.put({
            'error':self.return_code,
            'local_min_structure':self.local_min_struc,
            'last_structure':self.last_struc,
            'alpha':self.alpha,
            'index':index,
            'time':dt
        })

    def start_XTB(self,queue,index):
        DIR = os.path.dirname(self.INP)
        os.chdir(DIR)
        pat = re.compile('ENERGY[=|\s]*(\S+)\s*VBIAS[=|\s]*(\S*)\s*MD_STEP[=|\s]*(\S*)')
        err = open(self.ERR,'a')
        t0 = time.time()
        popen = subprocess.Popen(self.CMD_LINE, stdout=subprocess.PIPE, universal_newlines=True,stderr=err)
        for stdout_line in iter(popen.stdout.readline, ""):
            res = re.search(pat,stdout_line)
            if res != None:
                e,vb,n = re.findall(pat,stdout_line)[0]
                try:
                    e = float(e)
                except:
                    e = float('NaN')
                try:
                    vb = float(vb)
                except:
                    vb = float('NaN')
                try:
                    n = int(n)
                except:
                    n = float('NaN')
                self.MD_steps+=1
                self.ener_shift.append(e)
                self.vb_shift.append(vb)
        self.return_code = popen.wait()
        popen.stdout.close()
        err.close()
        t = time.time()
        if self.return_code:
            with open(self.INP,'r') as fi:
                err_inp = fi.read()
            with open(self.XYZ,'r') as fi:
                err_xyz = fi.read()
            err = open(self.ERR,'a')
            err.write(f'\ncalculation is down at {self.MD_steps} step\n')
            err.write(f'command line args: {self.CMD_LINE}\n')
            err.write(f'input file: \n{err_inp}\n')
            err.write(f'xyz file: \n{err_xyz}\n')
            err.close()
        else:
            assert len(self.vb_shift) == len(self.ener_shift) == self.MD_steps
            # processing local minimum structure
            min_e = self.ener_shift[0]
            min_i = 0
            for i,e in enumerate(self.ener_shift):
                if e < min_e:
                    min_e = e
                    min_i = i
            min_vb = self.vb_shift[min_i]
            min_xyz = base_utils.get_frame_xyz(self.TRJ,min_i)
            self.local_min_struc['Energy'] = min_e
            self.local_min_struc['Vbias'] = min_vb
            self.local_min_struc['xyz'] = min_xyz
            self.local_min_struc['index'] = min_i
            self.local_min_struc['Nsteps'] = self.MD_steps
            # processing last structure
            last_e = self.ener_shift[-1]
            last_vb = self.vb_shift[-1]
            last_i = self.MD_steps-1
            last_xyz = base_utils.get_frame_xyz(self.TRJ)
            self.last_struc['Energy'] = last_e
            self.last_struc['Vbias'] = last_vb
            self.last_struc['xyz'] = last_xyz
            self.last_struc['index'] = last_i
        os.chdir('..')
        dt = t-t0
        queue.put({
            'error':self.return_code,
            'local_min_structure':self.local_min_struc,
            'last_structure':self.last_struc,
            'alpha':self.alpha,
            'index':index,
            'time':dt
        })

    def start_LAMMPS(self,queue,index):
        header_found = False
        DIR = self.DIR
        os.chdir(DIR)
        pat = re.compile('Time\s*Temp\s*PotEng\s*KinEng\s*Press\s*Volume\s*Density')
        err = open(self.ERR,'a')
        t0 = time.time()
        popen = subprocess.Popen(self.CMD_LINE, stdout=subprocess.PIPE, universal_newlines=True,stderr=err)
        for stdout_line in iter(popen.stdout.readline, ""):
            res = re.search(pat,stdout_line)
            if res != None:
                header_found = True
                continue
            if re.search('Performance',stdout_line) !=None and header_found:
                break
            if header_found:
                try:
                    _, _, e, _, _, _, _ = [int(stdout_line.split()[0])] + [float(i) for i in stdout_line.split()[1:]]
                    vb = float('NaN')
                except:
                    continue
                self.MD_steps+=1
                self.ener_shift.append(e)
                self.vb_shift.append(vb)
        self.return_code = popen.wait()
        popen.stdout.close()
        err.close()
        t = time.time()
        if self.return_code:
            with open(self.INP,'r') as fi:
                err_inp = fi.read()
            # with open(self.XYZ,'r') as fi:
            #     err_xyz = fi.read()
            err = open(self.ERR,'a')
            err.write(f'\ncalculation is down at {self.MD_steps} step\n')
            err.write(f'command line args: {self.CMD_LINE}\n')
            err.write(f'input file: \n{err_inp}\n')
            # err.write(f'xyz file: \n{err_xyz}\n')
            err.close()
        else:
            assert len(self.vb_shift) == len(self.ener_shift) == self.MD_steps
            # processing local minimum structure
            min_e = self.ener_shift[0]
            min_i = 0
            for i,e in enumerate(self.ener_shift):
                if e < min_e:
                    min_e = e
                    min_i = i
            min_vb = self.vb_shift[min_i]
            min_xyz = base_utils.get_frame_xyz(self.TRJ,min_i)
            self.local_min_struc['Energy'] = min_e
            self.local_min_struc['Vbias'] = min_vb
            self.local_min_struc['xyz'] = min_xyz
            self.local_min_struc['index'] = min_i
            self.local_min_struc['Nsteps'] = self.MD_steps
            # processing last structure
            last_e = self.ener_shift[-1]
            last_vb = self.vb_shift[-1]
            last_i = self.MD_steps-1
            last_xyz = base_utils.get_frame_xyz(self.TRJ)
            self.last_struc['Energy'] = last_e
            self.last_struc['Vbias'] = last_vb
            self.last_struc['xyz'] = last_xyz
            self.last_struc['index'] = last_i
        os.chdir('..')
        dt = t-t0
        queue.put({
            'error':self.return_code,
            'local_min_structure':self.local_min_struc,
            'last_structure':self.last_struc,
            'alpha':self.alpha,
            'index':index,
            'time':dt
        })




# os.chdir('/home/artem/LAMMPS_TEST/dbg/')

# w1 = MD_calc('/home/artem/LAMMPS_TEST/dbg/',None,None,[os.path.join(LAMMPS_PATH,'lmp'),'-i','lammps.inp'],
#              'lammps.inp','OUT','ERR','traj.xyz','lammps.restart',False,
#              1.0,0,'structures.xyz','lammps')


# print(w1.local_min_struc)

# w1.start_LAMMPS2()
# print(w1.local_min_struc)

# print(0)

# TODO
# METHODS TO PERFORM TEST SITUATION: known structure is formed in World X at step Y

# 1. Separate trajectory of replica evolution during next steps

# 2. Probability of exchange statistics
        
# 3. Visualisation in Pyplot
        

# TEST 2: H-bond freezes




class World:
    WALL = False
    Ncycles = 0
    error = 0
    struc_before_exchange = None
    struc_after_exchange = None
    CMD_LINE = None
    exchanged = False
    shake = False
    external_pot = True


    def __init__(self, BASE_DIR, alpha, TOP, XYZ, INP, TRJ, RST, RSTBIN, XYZF, select_last_frame, engine) -> None:
        assert os.path.isfile(os.path.join(BASE_DIR,TOP)), 'parameter file not found!'
        assert os.path.isfile(os.path.join(BASE_DIR,INP)), 'input template not found!'
        assert engine in ('cp2k','xtb','lammps'), f'MD engine not found: {engine}!'
        try:
            alpha = float(alpha)
        except:
            raise TypeError(f'could not convert alpha to Float: {alpha}')
        self.alpha = alpha
        self.engine = engine
        self.BASE_DIR = BASE_DIR
        self.DIR = os.path.join(BASE_DIR,f'{alpha:.3e}')
        try:
            os.mkdir(self.DIR)
        except:
            shutil.rmtree(self.DIR, ignore_errors=True)
            os.mkdir(self.DIR)
        self.TOP = os.path.join(self.DIR,TOP)
        self.TRJ = os.path.join(self.DIR,TRJ)
        if RST != None:
            self.RST = os.path.join(self.DIR,RST)
            shutil.copyfile(os.path.join(BASE_DIR,RST), os.path.join(self.DIR,RST))
        else:
            self.RST  = None
        self.RSTBIN = os.path.join(self.DIR,RSTBIN)
        self.INP = INP
        self.OUT = INP[:-3]+'out'
        self.ERR = INP[:-3]+'err'
        self.XYZ_FOUND = os.path.join(self.BASE_DIR,XYZF)
        self.select_last_frame = select_last_frame
        if XYZ != None:
            self.XYZ = os.path.join(self.DIR,XYZ)
            shutil.copyfile(os.path.join(BASE_DIR,XYZ), os.path.join(self.DIR,self.XYZ))
        else:
            self.XYZ = None
        shutil.copyfile(os.path.join(BASE_DIR,TOP), os.path.join(self.DIR,self.TOP))
        shutil.copyfile(os.path.join(BASE_DIR,INP), os.path.join(self.DIR,self.INP))
        self.ARGS = None
        if self.engine == 'cp2k.inp':
            self.ARGS = base_utils.read_cp2k()
            self.ARGS['FORCE_EVAL'][0]['MM'][0]['FORCEFIELD'][0]['FORCE_SCALE'][0] = str(self.alpha)
            self.CMD_LINE = [os.path.join('cp2k.sopt'),'-i',self.INP]
        elif self.engine == 'xtb':
            self.CMD_LINE = [os.path.join(XTB_PATH,'xtb'),self.XYZ,'-I',self.INP,'--md','--hremd','--gfnff','-P','1']
            if os.path.exists(self.XYZ_FOUND) and '--xyzset' not in self.CMD_LINE:
                self.CMD_LINE.append('--xyzset')
                self.CMD_LINE.append(self.XYZ_FOUND)
            if not self.WALL:
                self.ARGS.pop('$wall','')
        elif self.engine == 'lammps':
            self.CMD_LINE = [os.path.join(LAMMPS_PATH,'lmp'),'-i',self.INP]

    def remove_ext_pot(self):
        if self.engine != 'cp2k':
            print(f'External potential was not removed in a={self.alpha}: not implemented for {self.engine} engine!')
            return
        self.ARGS['FORCE_EVAL'][0].pop('EXTERNAL_POTENTIAL','')
        self.external_pot = False
        print(f'External potential has been removed in a={self.alpha}')

    def set_pot(self):
        if self.engine != 'cp2k':
            print(f'External potential was not set in a={self.alpha}: not implemented for {self.engine} engine!')
            return
        # results of MDs with different parameter sets (a,b,c):
        # 1. 1000,50,5: too weak
        # 1. 2000,60,10: too weak
        # 1. 4000,70,20: too strong

        N=10
        C = 10**-N
        a=1000
        b=50
        c=10
        fx = f"{C:.{N}f}*(X^2)^4 - {a}*exp(-((X-{b})^2)/{c}^2)"
        fy = f"{C:.{N}f}*(Y^2)^4 - {a}*exp(-((Y-{b})^2)/{c}^2)"
        fz = f"{C:.{N}f}*(Z^2)^4 - {a}*exp(-((Z-{b})^2)/{c}^2)"
        self.ARGS['FORCE_EVAL'][0]['EXTERNAL_POTENTIAL'][0]['FUNCTION'][0] = fz
        self.ARGS['FORCE_EVAL'][0]['EXTERNAL_POTENTIAL'][1]['FUNCTION'][0] = fy
        self.ARGS['FORCE_EVAL'][0]['EXTERNAL_POTENTIAL'][2]['FUNCTION'][0] = fx

# move center of external potential
    def shake_cell(self,N,M=10,dr=20):
        if self.engine != 'cp2k':
            print(f'No shaking cell a={self.alpha}: not implemented for {self.engine} engine!')
            return
        if not self.external_pot:
            return
        C = 10**-M
        if N%12==1:
            fz = f"{C:.{M}f}*((Z-{dr})^2)^4"
            fy = f"{C:.{M}f}*(Y^2)^4"
            fx = f"{C:.{M}f}*(X^2)^4"
        elif N%12==3:
            fz = f"{C:.{M}f}*((Z+{dr})^2)^4"
            fy = f"{C:.{M}f}*(Y^2)^4"
            fx = f"{C:.{M}f}*(X^2)^4"
        elif N%12==5:
            fz = f"{C:.{M}f}*(Z^2)^4"
            fy = f"{C:.{M}f}*((Y-{dr})^2)^4"
            fx = f"{C:.{M}f}*(X^2)^4"
        elif N%12==7:
            fz = f"{C:.{M}f}*(Z^2)^4"
            fy = f"{C:.{M}f}*((Y+{dr})^2)^4"
            fx = f"{C:.{M}f}*(X^2)^4"
        elif N%12==9:
            fz = f"{C:.{M}f}*(Z^2)^4"
            fy = f"{C:.{M}f}*(Y^2)^4"
            fx = f"{C:.{M}f}*((X-{dr})^2)^4"
        elif N%12==11:
            fz = f"{C:.{M}f}*(Z^2)^4"
            fy = f"{C:.{M}f}*(Y^2)^4"
            fx = f"{C:.{M}f}*((X+{dr})^2)^4"
        else:
            fz = f"{C:.{M}f}*(Z^2)^4"
            fy = f"{C:.{M}f}*(Y^2)^4"
            fx = f"{C:.{M}f}*(X^2)^4"
        self.ARGS['FORCE_EVAL'][0]['EXTERNAL_POTENTIAL'][0]['FUNCTION'][0] = fz
        self.ARGS['FORCE_EVAL'][0]['EXTERNAL_POTENTIAL'][1]['FUNCTION'][0] = fy
        self.ARGS['FORCE_EVAL'][0]['EXTERNAL_POTENTIAL'][2]['FUNCTION'][0] = fx
    
    def set_results(self,res):
        self.error = res['error']
        assert self.alpha == res['alpha'], f'Error in Queue order: {res["alpha"]} --> {self.alpha}'
        if self.error:
            # TODO
            self.struc_before_exchange = {
                'Energy':float('NaN'),
                'Vbias':float('NaN'),
                'xyz':None,
                'index':None,
                'Nsteps':None}
        else:
            # skip structure with min.E if it appeared too early
            skip_local = (res['local_min_structure']['index'] < res['local_min_structure']['Nsteps']/2)
            if self.select_last_frame | skip_local:
                self.struc_before_exchange = res['last_structure'].copy()
            else:
                self.struc_before_exchange = res['local_min_structure'].copy()
        self.Ncycles+=1
            
    def get_structure_before(self):
        return copy.deepcopy(self.struc_before_exchange)
    
    def get_structure_after(self):
        return copy.deepcopy(self.struc_after_exchange)
    
    def set_structure_after(self,struc):
        self.struc_after_exchange = struc

    def update_xyz(self):
        if self.engine != 'lammps':
            with open(self.XYZ,'w') as xyz:
                xyz.write(self.struc_after_exchange['xyz'])
        else:
            raise RuntimeError(f'not implemented for LAMMPS!')

    def get_binrestart_name(self):
        try:
            return os.path.abspath(self.RSTBIN)
        except:
            raise RuntimeError(f'Restart files are not implemented for "{self.engine}"!')

    def update_shared_xyz(self):
        if os.path.exists(self.XYZ_FOUND):
            with open(self.XYZ_FOUND,'a') as xyz:
                xyz.write('\n'+self.struc_after_exchange['xyz'])
        else:
            with open(self.XYZ_FOUND,'w') as xyz:
                xyz.write(self.struc_after_exchange['xyz'])

        # if '--xyzset' not in self.CMD_LINE:
        #     self.CMD_LINE.append('--xyzset')
        #     self.CMD_LINE.append(self.XYZ_FOUND)

    def rename_trj_file(self):
        dirname = os.path.dirname(self.TRJ)
        fname = os.path.basename(self.TRJ)
        assert fname.count('.') == 1, 'change condition in this method'
        base,ext = fname.split('.')
        fname_new = f'{base}-{self.Ncycles}.{ext}'
        new = os.path.join(dirname,fname_new)
        try:
            shutil.move(self.TRJ, new)
        except:
            #TODO: possible case if calculation is down
            raise RuntimeError(f'could not find trajectory file: {self.TRJ}')

