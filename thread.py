import re,os,subprocess,shutil,time,copy
import xtb_IO
import results
# USE_RMSD = False
# RMSD_REF = '/home/artem/HREMD_algorithm/cp2k/rmsd/second.xyz'
# RMSD_REF = '/home/artem/HREMD_algorithm/cp2k/rmsd/atropo2-cp2kopt.xyz'
# REF_XYZ = results.get_structure_xyz(RMSD_REF)


CP2K = True

XTB_PATH = '/home/artem/xtb-v2/build'
BLAS_PATH = '/home/artem/OpenBLAS/lib/'
# CP2K_PATH = '/home/artem/cp2k/exe/local/'
CP2K_PATH = ''
# TRJ_NAME = 'xtb.trj'
TRJ_NAME = 'PAM-pos-1.xyz'
# TRJ_NAME = 'Cu2L2_coord-pos-1.xyz'
# INP_NAME = 'xtb.inp'
INP_NAME = 'cp2k.inp'
# ERR_NAME = 'xtb.err'
ERR_NAME = 'cp2k.err'
# OUT_NAME = 'xtb.out'
OUT_NAME = 'cp2k.out'
# single MD calculation; lives while xtb runs


GEO_OPT_PARM = {
    "GLOBAL": [
        {
            "PRINT_LEVEL": [
                "SILENT"
            ],
            "PROJECT": [
            f"CP2K_opt"
            ],
            "RUN_TYPE": [
                "GEO_OPT"
            ]
        }
    ],
    "FORCE_EVAL": [
        {
            "METHOD": [
                "FIST"
            ],
            "MM": [
                {
                    "FORCEFIELD": [
                        {
                            "DO_NONBONDED": [
                                "T"
                            ],
                            "ALL_EI_SCALE": [
                                "1.0"
                            ],
                            "ALL_VDW_SCALE": [
                                "1.0"
                            ],
                            "EI_SCALE14": [
                                "1.0"
                            ],
                            "VDW_SCALE14": [
                                "1.0"
                            ],
                            "FORCE_SCALE": [
                                "1.0"
                            ],
                            "PARMTYPE": [
                                "AMBER"
                            ],
                            "PARM_FILE_NAME": [
                                ""
                            ],
                            "SHIFT_CUTOFF": [
                                "F"
                            ],
                            "SPLINE": [
                                {
                                    "EMAX_SPLINE": [
                                        "1.0E+08"
                                    ],
                                    "RCUT_NB": [
                                        "[angstrom] 1.0E+01"
                                    ]
                                }
                            ]
                        }
                    ],
                    "POISSON": [
                        {
                            "EWALD": [
                                {
                                    "ALPHA": [
                                        ".40"
                                    ],
                                    "EWALD_TYPE": [
                                        "NONE"
                                    ],
                                    "GMAX": [
                                        "80"
                                    ]
                                }
                            ],
                            "PERIODIC": [
                                "NONE"
                            ]
                        }
                    ]
                }
            ],
            "RESCALE_FORCES": [
                {
                    "MAX_FORCE": [
                        "10"
                    ]
                }
            ],
            "STRESS_TENSOR": [
                "NONE"
            ],
            "SUBSYS": [
                {
                    "CELL": [
                        {
                            "ABC": [
                                "[angstrom] 100 100 100"
                            ],
                            "ALPHA_BETA_GAMMA": [
                                "90 90 90"
                            ],
                            "PERIODIC": [
                                "NONE"
                            ]
                        }
                    ],
                    "TOPOLOGY": [
                        {
                            "CENTER_COORDINATES": [
                                {
                                    "CENTER_POINT": [
                                        "0 0 0"
                                    ]
                                }
                            ],
                            "CONN_FILE_FORMAT": [
                                "AMBER"
                            ],
                            "CONN_FILE_NAME": [
                                ""
                            ],
                            "COORD_FILE_FORMAT": [
                                "XYZ"
                            ],
                            "COORD_FILE_NAME": [
                                ""
                            ]
                        }
                    ]
                }
            ]
        }
    ],
    "MOTION": [
        {
            "GEO_OPT": [
                {
                    "MAX_ITER": [
                        "5000"
                    ]
                }
            ],
            "PRINT": [
                {
                    "RESTART": [
                        {
                            "BACKUP_COPIES": [
                                "0"
                            ],
                            "EACH": [
                                {
                                    "GEO_OPT": [
                                        "0"
                                    ]
                                }
                            ]
                        }
                    ],
                    "RESTART_HISTORY": [
                        {
                            "EACH": [
                                {
                                    "GEO_OPT": [
                                        "0"
                                    ]
                                }
                            ]
                        }
                    ],
                    "TRAJECTORY": [
                        {
                            "EACH": [
                                {
                                    "MD": [
                                        "1"
                                    ]
                                }
                            ],
                            "FORMAT": [
                                "XYZ"
                            ]
                        }
                    ]
                }
            ]
        }
    ]
}

def start_CP2K(queue,index,INP,ERR,XYZ,TRJ,CMD_LINE,alpha):
    MD_steps = 0
    ener_shift = []
    vb_shift = []

    local_min_struc = {
        'Energy':None,
        'Vbias':None,
        'xyz':None
        }
    last_struc = {
        'Energy':None,
        'Vbias':None,
        'xyz':None
        }

    DIR = os.path.dirname(INP)
    os.chdir(DIR)
    pat = re.compile('ENERGY[=|\s]*(\S+)\s*VBIAS[=|\s]*(\S*)\s*MD_STEP[=|\s]*(\S*)')
    err = open(ERR,'a')
    t0 = time.time()
    popen = subprocess.Popen(CMD_LINE, stdout=subprocess.PIPE, universal_newlines=True,stderr=err)
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
            MD_steps+=1
            ener_shift.append(e)
            vb_shift.append(vb)
    return_code = popen.wait()
    popen.stdout.close()
    err.close()
    t = time.time()
    if return_code:
        with open(INP,'r') as fi:
            err_inp = fi.read()
        with open(XYZ,'r') as fi:
            err_xyz = fi.read()
        err = open(ERR,'a')
        err.write(f'\ncalculation is down at {MD_steps} step\n')
        err.write(f'command line args: {CMD_LINE}\n')
        err.write(f'input file: \n{err_inp}\n')
        err.write(f'xyz file: \n{err_xyz}\n')
        err.close()
    else:
        assert len(vb_shift) == len(ener_shift) == MD_steps
        # processing local minimum structure
        min_e = ener_shift[0]
        min_i = 0
        for i,e in enumerate(ener_shift):
            if e < min_e:
                min_e = e
                min_i = i
        min_vb = vb_shift[min_i]
        min_xyz = xtb_IO.get_frame_xyz(TRJ,min_i)
        local_min_struc['Energy'] = min_e
        local_min_struc['Vbias'] = min_vb
        local_min_struc['xyz'] = min_xyz
        local_min_struc['index'] = min_i
        local_min_struc['Nsteps'] = MD_steps
        # processing last structure
        last_e = ener_shift[-1]
        last_vb = vb_shift[-1]
        last_i = MD_steps-1
        last_xyz = xtb_IO.get_frame_xyz(TRJ)
        last_struc['Energy'] = last_e
        last_struc['Vbias'] = last_vb
        last_struc['xyz'] = last_xyz
        last_struc['index'] = last_i
    os.chdir('..')
    dt = t-t0
    queue.put({
        'error':return_code,
        'local_min_structure':local_min_struc,
        'last_structure':last_struc,
        'alpha':alpha,
        'index':index,
        'time':dt
    })


class MD_calc:
    local_min_struc = {
        'Energy':None,
        'Vbias':None,
        'xyz':None
        }
    last_struc = {
        'Energy':None,
        'Vbias':None,
        'xyz':None
        }
    min_rmsd_struc = {
        'RMSD':None,
        'xyz':None,
        'Energy':None,
        'Vbias':None
        }
    ener_shift = []
    rmsd_shift = []
    vb_shift = []
    MD_steps = 0
    error_code = 0

    def __init__(self,XYZ,ARGS,CMD_LINE,INP,OUT,ERR,TRJ,alpha,Nxyz,XYZ_SHARED) -> None:

        # try:
        #     lib_omp = os.environ['OPENMPI_LIBS']
        # except:
        #     raise EnvironmentError('type in console: source /home/artem/cp2k/tools/toolchain/install/setup')
        if not(os.path.exists(ERR)):
            open(ERR,'w').close()
        if not(os.path.exists(OUT)):
            open(OUT,'w').close()
        # parm = xtb_IO.update_parm(ARGS)
        # # TODO: move this part inside of xtb
        # if Nxyz != 0:
        #     parm['$metadyn']['save'] = Nxyz
        #     if '--xyzset' not in CMD_LINE:
        #         assert os.path.exists(XYZ_SHARED)
        #         CMD_LINE.append('--xyzset')
        #         CMD_LINE.append(XYZ_SHARED)
        # parm.pop('$wall','')
        # xtb_IO.write_input(parm,INP)
        xtb_IO.write_input_cp2k(ARGS,fname=INP)
        self.INP = INP
        self.OUT = OUT
        self.ERR = ERR
        self.XYZ = XYZ
        self.TRJ = TRJ
        self.ARGS = ARGS # dict with sections of input file
        self.CMD_LINE = CMD_LINE # [path/to/xtb, '--arg1', '--arg2', ...]
        self.alpha = alpha





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
            min_xyz = xtb_IO.get_frame_xyz(self.TRJ,min_i)
            self.local_min_struc['Energy'] = min_e
            self.local_min_struc['Vbias'] = min_vb
            self.local_min_struc['xyz'] = min_xyz
            self.local_min_struc['index'] = min_i
            self.local_min_struc['Nsteps'] = self.MD_steps
            # processing last structure
            last_e = self.ener_shift[-1]
            last_vb = self.vb_shift[-1]
            last_i = self.MD_steps-1
            last_xyz = xtb_IO.get_frame_xyz(self.TRJ)
            self.last_struc['Energy'] = last_e
            self.last_struc['Vbias'] = last_vb
            self.last_struc['xyz'] = last_xyz
            self.last_struc['index'] = last_i
            # if USE_RMSD:
            #     # self.rmsd_shift = results.rmsd_xyz(self.TRJ, REF_XYZ)
            #     self.rmsd_shift = results.rmsd_pyxyz(self.TRJ, RMSD_REF)
            #     min_rmsd = self.rmsd_shift[0]
            #     min_i = 0
            #     for i,val in enumerate(self.rmsd_shift):
            #         if val < min_rmsd:
            #             min_rmsd = val
            #             min_i = i
            #     if (min_i < self.MD_steps/5):
            #         min_i = self.MD_steps-1
            #         min_rmsd_xyz = xtb_IO.get_frame_xyz(self.TRJ,min_i)
            #     else:
            #         min_rmsd_xyz = xtb_IO.get_frame_xyz(self.TRJ,min_i)
            #     min_vb = self.vb_shift[min_i]
            #     min_e = self.ener_shift[min_i]
            #     self.min_rmsd_struc['Energy'] = min_e
            #     self.min_rmsd_struc['Vbias'] = min_vb
            #     self.min_rmsd_struc['RMSD'] = min_rmsd
            #     self.min_rmsd_struc['xyz'] = min_rmsd_xyz
            #     self.min_rmsd_struc['Nsteps'] = self.MD_steps
            #     self.min_rmsd_struc['index'] = min_i
        os.chdir('..')
        dt = t-t0
        queue.put({
            'error':self.return_code,
            'local_min_structure':self.local_min_struc,
            'rmsd_structure':self.min_rmsd_struc,
            'last_structure':self.last_struc,
            'alpha':self.alpha,
            'index':index,
            'time':dt
        })


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


    def __init__(self,BASE_DIR,alpha,T_REMD,TOP,XYZ,XYZ_FOUND,ARGS,USE_LAST_STRUC) -> None:
        assert os.path.isfile(os.path.join(BASE_DIR,TOP)), '.prmtop not found!'
        try:
            alpha = float(alpha)
        except:
            raise TypeError(f'could not convert alpha to Float: {alpha}')
        self.alpha = alpha
        self.T_REMD = T_REMD
        self.BASE_DIR = BASE_DIR
        self.DIR = os.path.join(BASE_DIR,'{:.3e}'.format(alpha))
        try:
            os.mkdir(self.DIR)
        except:
            shutil.rmtree(self.DIR, ignore_errors=True)
            os.mkdir(self.DIR)
        self.TOP = os.path.join(self.DIR,TOP)
        self.XYZ = os.path.join(self.DIR,XYZ)
        self.TOP = os.path.join(self.DIR,TOP)
        self.TRJ = os.path.join(self.DIR,TRJ_NAME)
        self.INP = os.path.join(self.DIR,INP_NAME)
        self.OUT = os.path.join(self.DIR,OUT_NAME)
        self.ERR = os.path.join(self.DIR,ERR_NAME)
        self.XYZ_FOUND = os.path.join(self.BASE_DIR,XYZ_FOUND)
        shutil.copyfile(os.path.join(BASE_DIR,TOP), os.path.join(self.DIR,self.TOP))
        # shutil.copyfile(os.path.join(BASE_DIR,TOP), os.path.join(self.DIR,'reference.xyz'))
        shutil.copyfile(os.path.join(BASE_DIR,XYZ), os.path.join(self.DIR,self.XYZ))
        self.ARGS = copy.deepcopy(ARGS) # dict
        # self.ARGS['FORCE_EVAL']['MM']['FORCEFIELD']['FORCE_SCALE'] = str(self.alpha)
        if self.T_REMD:
            self.ARGS['FORCE_EVAL'][0]['MM'][0]['FORCEFIELD'][0]['FORCE_SCALE'][0] = str(1.0)
            self.ARGS['MOTION'][0]['MD'][0]['TEMPERATURE'][0] = str(self.alpha)
        else:
            self.ARGS['FORCE_EVAL'][0]['MM'][0]['FORCEFIELD'][0]['FORCE_SCALE'][0] = str(self.alpha)
        self.USE_LAST_STRUC = USE_LAST_STRUC
        # if int(self.alpha) == 1:
            # self.ARGS.pop('$constrain','')
        if not self.WALL:
            self.ARGS.pop('$wall','')
        # self.CMD_LINE = [os.path.join(XTB_PATH,'xtb'),self.XYZ,'-I',self.INP,'--md','--hremd','--gfnff','-P','1','--norestart']
        # self.CMD_LINE = [os.path.join(XTB_PATH,'xtb'),self.XYZ,'-I',self.INP,'--md','--hremd','--gfnff','-P','1']
        self.CMD_LINE = [os.path.join(CP2K_PATH,'cp2k.sopt'),'-i',self.INP]
        # if os.path.exists(self.XYZ_FOUND) and '--xyzset' not in self.CMD_LINE:
            # self.CMD_LINE.append('--xyzset')
            # self.CMD_LINE.append(self.XYZ_FOUND)

    def remove_ext_pot(self):
        self.ARGS['FORCE_EVAL'][0].pop('EXTERNAL_POTENTIAL','')
        self.external_pot = False
        print(f'External potential has been removed in a={self.alpha}')

    def set_pot(self):
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




    def start_MD(self,queue,i):
        # print(self.alpha)
        calc = MD_calc(self.XYZ,self.ARGS,self.CMD_LINE,self.INP,self.OUT,self.ERR,self.TRJ,self.alpha,i)
        res = calc.start_XTB()
        queue.put(res)
        # print(res['error'])
    
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
            # if USE_RMSD:
            #     self.struc_before_exchange = res['rmsd_structure'].copy()
            #     self.Ncycles+=1
            #     return
            skip_local = (res['local_min_structure']['index'] < res['local_min_structure']['Nsteps']/2)
            if self.USE_LAST_STRUC | skip_local:
                self.struc_before_exchange = res['last_structure'].copy()
                # print(f'LOCAL MIN ATTEMPT: {res["local_min_structure"]["index"]} - DISCARDED')
            else:
                # print(f'LOCAL MIN ATTEMPT: {res["local_min_structure"]["index"]} - ACCEPTED')
                self.struc_before_exchange = res['local_min_structure'].copy()
        self.Ncycles+=1
            
    def get_structure_before(self):
        return copy.deepcopy(self.struc_before_exchange)
    
    def get_structure_after(self):
        return copy.deepcopy(self.struc_after_exchange)
    
    def set_structure_after(self,struc):
        self.struc_after_exchange = struc

    def update_xyz(self):
        with open(self.XYZ,'w') as xyz:
            xyz.write(self.struc_after_exchange['xyz'])
    
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

