import os,shutil, time, tempfile
import numpy as np
import multiprocessing as mp
import base_utils
import engines

# TODO добавить возможность обновлять координаты через read_dump
# TODO рестарт упавшего/завершенного расчета через сериализацию
# TODO отладить удаление миров

''' Representation of crystal structure'''
class Structure:
    def __init__(self,topology, a, b, c, alpha, beta, gamma, xyz=None, symm = 'P1') -> None:
        self.e = None

        self.topology = topology
        # self.xyz = xyz
        if xyz == None:
            self.xyz = base_utils.getXyzLmp(topology)
            # self.xyz = base_utils.getXyzfromData(topology)
        else:
            self.xyz = xyz
        self.a = a
        self.b = b
        self.c = c
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.symm = symm

    def _calcPxrd(self):
        pass

    def calcPxrd(self,ref):
        pass

    def calcRmsd(self,ref):
        pass

    def getXyz(self,internal=False):
        assert self.xyz != None
        return self.xyz
    
    def changeXyz(self,xyz, energy=None, internal=False):
        assert type(xyz) is list
        assert len(xyz) == len(self.xyz), 'number of atoms does not match' # TODO: check order of atoms
        # assert len(xyz.split('\n')) == len(self.xyz.split('\n')), 'number of atoms does not match' # TODO: check order of atoms
        self.xyz = xyz
        self.e = energy

    def changeCell(self):
        pass

    def getEnergy(self):
        if self.e == None:
            return float('Nan')
        return self.e

    def relax(self, engine, method, maxiter): # TODO
        pass

    def _updateTopo(self):
        pass

    def equilibrate(self, WD, datafile=None,parm=None, vars=None):
        return 'Not implemented'

    def getContacts(self): # с этим несколько сложнее
        pass

    def guessSym(self): # еще сложнее
        pass

    def reduceSpaceGroup(self):
        pass

    def changeSpaceGroup(self):
        pass

    def setConstraints(self):
        pass

    def getConstraints(self):
        pass

'''Object for collecting statistics on exchanges'''
class ExchangeCounter:

    def __init__(self, alphaSet, WD):
        self.indices = [i for i in range(len(alphaSet))]
        self.imap = {i:self.indices.index(i) for i in self.indices}
        self.rows = []
        self.deleted = []
        self.dumpFname = os.path.join(WD,'exchanges.log')
        
        open(self.dumpFname,'w').close()

    def addIndex(self, j):
        # добавление происходит справа: т.е. i-й элемент < j остается при своем значении, 
        # j становится на j-е место, а элементы после j смещаются на +1; при этом первый
        # и последний элементы всегда остаются на своих местах
        assert 0 < j <= len(self.indices)-1, 'wrong index'
        for i in range(len(self.indices)):
            if self.indices[i]>=j:
                self.indices[i]+=1
        for row in self.rows:
            for i in range(len(self.indices)):
                if row[i]>=j:
                    row[i]+=1
            row.insert(j,float('Nan'))
        self.indices.insert(j,j)
        self.imap = {i:i for i in range(len(self.indices))}

    def delIndex(self, j): # TODO: debug this
        assert 0 < j < len(self.indices)-1, 'wrong index'
        val = self.indices.pop(j)
        for i in range(len(self.indices)):
            if self.indices[i]>val:
                self.indices[i]-=1
        self.imap = {i:i for i in range(len(self.indices))}
        self.deleted.append(j)

    def updateIndex(self, i, j):
        assert i<j, 'unexpected order of indexes'
        i1 = self.indices.index(i)
        i2 = self.indices.index(j)
        atmp = self.indices[i1]
        self.indices[i1] = self.indices[i2]
        self.indices[i2] = atmp

    def addRow(self):
        row = [self.imap[i] for i in self.indices]
        for i in reversed(self.deleted):
            for j in range(len(row)):
                if row[j]>=i:
                    row[j]+=1
            row.insert(i, float('Nan'))
        self.rows.append(row)
    
    def writeLog(self):
        with open(self.dumpFname,'a') as fo:
            for row in self.rows:
                fo.write('   '.join([str(i) for i in row]) + '\n')


'''
Representation of REMD simulation

 - `alphaSet`: list : [ (alpha(1), Structure(1), Simulation(1), *other), (alpha(2), Structure(2), Simulation(2), *other) ... ];
 - `swapCount` : list : [n1, n2, ...];
 - `Counter` : ExchangeCounter : ExchangeCounter(basedir)
 - `initStruc` : Structure : Structure(datafile, a, b, c, alpha, beta, gamma)

'''
class REMD:
    unrestrictedExchange = True
    changeOrder = False
    # selectLastStruc = False
    addWorlds = True
    delWorlds = False
    Ncheck = 10 # number of iterations between attempts to change alpha set
    Pmin = 10 # в %
    Pmax = 80 # в %

    

    # assert selectLastStruc == True, 'search for local minimum is not implemented'
    assert delWorlds is False, 'deleting worlds is not implemented'

    def __init__(self, alphaRange, datafile, baseDir, Nmax=15, seed=999999, T=273.15, NPTs=[], delayNPT=1, Ncores=1, selectLastStruc=True,  parm=None, vars=None, options=None, *other):
        assert os.path.exists(datafile), f'Data file not found: {datafile}'
        assert len(alphaRange) > 1, 'Not enough initial worlds for REMD (N must be >=2)'
        assert Nmax >= len(alphaRange), f'The requested number of worlds ({len(alphaRange)}) exceeds Nmax = {Nmax}'
        try:
            os.makedirs(baseDir)
        except FileExistsError:
            print(f'Directory already exists: {baseDir}')
            shutil.rmtree(baseDir)
            os.makedirs(baseDir)
            # raise # ради сохранности предыдуших расчетов
        self.baseDir = baseDir
        self.dumpXYZ = os.path.join(baseDir,'structures.xyz')
        self.dumpEnergy = os.path.join(baseDir,'energy.log')
        open(self.dumpXYZ,'w').close()
        open(self.dumpEnergy,'w').close()
        self.initDataFile = datafile
        self.initParm = parm
        self.initVars = vars
        self.initOptions = options
        self.initOther = other
        self.Nmax = Nmax
        self.seed = seed
        self.T = T
        self.selectLastStruc = selectLastStruc
        self.beta = 1/(T * 1.987204259 * 10**-3) # assert kcal/mol TODO: pass as global argument
        self.Nadd = 0
        self.Ndel = 0
        self.Niter = 0
        self.initStruc = Structure(datafile, None, None, None, None, None, None)
        self.alphaSet = self._buildAlphaSet(alphaRange, self.initStruc, datafile, parm=self.initParm, vars=self.initVars, options=self.initOptions, *other)
        self.idump = [i for i in range(len(self.alphaSet))]
        self.withNPT = [self.alphaSet[i][2] for i in NPTs]
        self.delayNPT = delayNPT
        
        # if NPTs != None: # TODO: more convinient and general way to do this?
        #     if type(NPTs) is int:
        #         assert NPTs < len(self.alphaSet), 'incorrect index'
        #         self.alphaSet[NPTs][2]._updateOpt({'doNPT': True})
        #         print(f'Flexible cell in world with alpha = {self.alphaSet[NPTs][0]}')
        #     elif type(NPTs) is list:
        #         for i in NPTs:
        #             self.alphaSet[i][2]._updateOpt({'doNPT': True})
        #             print(f'Flexible cell in world with alpha = {self.alphaSet[i][0]}')
        #     else:
        #         raise RuntimeError('unknown type, must be `list` or `integer`')
        self.Ncores = Ncores
        self.swapCount = [0 for _ in range(len(self.alphaSet)-1)]
        self.Counter = ExchangeCounter(self.alphaSet, baseDir)
        self.randomGenerator = np.random.default_rng(self.seed)

    def reset(self):
        # _buildAlphaSet
        # Counter
        # random
        # swapCount
        # + clear directories and log files
        pass


    def _buildAlphaSet(self, alphas, initStruc, datafile, parm=None, vars=None, options=None, *other):
        alphas = sorted(alphas)
        aSet = []
        wds = []
        nrep = 0
        for alpha in alphas:
            wd = os.path.join(self.baseDir, f'{alpha:.3e}')
            if wd in wds: # for the case with several identical alphas
                nrep+=1
                wd+=f'{nrep}'
            else:
                nrep=0
            wds.append(wd)
            sim = engines.Simulation(alpha, wd, datafile, parm=parm, vars=vars, options=options)
            s0 = initStruc
            xyz0 = s0.getXyz()
            struc = Structure(datafile, s0.a, s0.b, s0.c, s0.alpha, s0.beta, s0.gamma, xyz=xyz0) # TODO: make it more obvious
            world = (alpha, struc, sim, *other)
            aSet.append(world)
        print(f'Initial alpha Set was created from 1 parent structure {initStruc} and shared topology: {datafile}')
        return aSet

    def _runMDSingle(self, sim): # TODO: remove
        return sim._runLAMMPS()

    def _runMD(self):
        failed = []
        good = []
        sims = [i[2] for i in self.alphaSet] # not very careful TODO
        strucs = [i[1] for i in self.alphaSet] # not very careful TODO
        args = []
        for sim in sims:
            if sim in self.withNPT and self.Niter > self.delayNPT:
                print(f'Iter {self.Niter}: NPT enabled')
                sim.prepare(options={'doNPT':True})
            else:
                sim.prepare()
            args.append((sim.WD, sim.IN_FNAME, sim.vars['d'], sim.ERR_FNAME, sim.TRJ_FNAME, sim.k, sim.patt, self.Ncores))
        with mp.Pool(len(self.alphaSet)) as pool:
            # results = pool.map(self._runMDSingle, sims)
            results = pool.map(engines.runMD, args)
        for i, sim in enumerate(sims):
            sim.update(results[i])
        for i, struc in enumerate(strucs):
            if results[i]['exitCode']:
                if i==len(self.alphaSet)-1:
                    print('WARNING: LAMMPS failed in unbiased world')
                    # raise RuntimeError('ERROR termination: LAMMPS failed in unbiased world')
                failed.append(self.alphaSet[i])
            else:
                if self.selectLastStruc:
                    xyz = results[i]['lastStruc']['xyz']
                    e = results[i]['lastStruc']['energy']
                else:
                    xyz = results[i]['minEStruc']['xyz']
                    e = results[i]['minEStruc']['energy']

                struc.changeXyz(xyz,e)
                good.append(self.alphaSet[i])
        if len(failed) == len(self.alphaSet):
            raise RuntimeError('ERROR termination : All MDs failed')
        alpha1, _, sim1, *_ = good[-1]
        rst1 = sim1.getRestartFile()
        for w in failed:
            alpha, _, sim, *_ = w
            sim.updateRestartFile(rst1)
            print(f'Failed MD with alpha {alpha} will be restarted from World with alpha = {alpha1}')

    @staticmethod
    def calc_single_point(world, struc):
        with tempfile.TemporaryDirectory() as tmpdir:
            alpha, _, refsim, *_ = world
            # print(refsim.keys.copy())
            sim = engines.Simulation(alpha, tmpdir, refsim.datfile, options={'single_point':True}, parm=refsim.keys.copy())
            _xyz = struc.getXyz()
            if _xyz:
                sim.updateXyz(_xyz)
            sim.prepare()
            args = (sim.WD, sim.IN_FNAME, sim.vars['d'], sim.ERR_FNAME, sim.TRJ_FNAME, sim.k, sim.patt, 1)
            total_pe = engines.run_single_point(args)
            # sim.update()
            return total_pe

        
    def _calcP(self, w1, w2):
        a1, struc1, *_ = w1
        a2, struc2, *_ = w2
        classic_scheme = True # delta = beta*[-E(a1,q1) - E(a2,q2) + E(a2,q1) + E(a1,q2)]
        if classic_scheme:
            a1q1 = struc1.getEnergy()
            a2q2 = struc2.getEnergy()
            a1q2 = REMD.calc_single_point(w1,struc2)
            a2q1 = REMD.calc_single_point(w2,struc1)
            dE = a2q1+a1q2-a1q1-a2q2
            dbg=True
            if dbg:
                log = os.path.join(self.baseDir,'prob.txt')
                if not os.path.exists(log):
                    open(log,'w').close()
                _dE = - a1q1*a1 - a2q2*a2 + a2q2*a1 + a1q1*a2
                with open(log,'a') as fo:
                    fo.write(f'Iteration {self.Niter}: alpha1 = {a1}; alpha2 = {a2}\n')
                    fo.write(f'              delta = beta*[-E(a1,q1) - E(a2,q2) + E(a2,q1) + E(a1,q2)]\n')
                    fo.write(f'              delta = {self.beta}*[-{a1q1} - {a2q2} + {a2q1} + {a1q2}] = {dE*self.beta}\n')
                    fo.write(f'              p = {min(1.0, np.exp(-dE*self.beta))}\n')
                    fo.write(f'              _delta = beta*[-E(q1)*a1 - E(q2)*a2 + E(q2)*a1 + E(q1)*a2]\n')
                    fo.write(f'              _delta = {self.beta}*[-{a1q1*a1} - {a2q2*a2} + {a2q2*a1} + {a1q1*a2}] = {_dE*self.beta}\n')
                    fo.write(f'              _p = {min(1.0, np.exp(-_dE*self.beta))}\n')
                print(f'{self.Niter}: alpha1 = {a1}; alpha2 = {a2}')
                print(f'              delta = beta*[-E(a1,q1) - E(a2,q2) + E(a2,q1) + E(a1,q2)]')
                print(f'              delta = {self.beta}*[-{a1q1} - {a2q2} + {a2q1} + {a1q2}] = {dE*self.beta}')
                print(f'              p = {min(1.0, np.exp(-dE*self.beta))}')
                print(f'              _delta = beta*[-E(q1)*a1 - E(q2)*a2 + E(q2)*a1 + E(q1)*a2]')
                print(f'              _delta = {self.beta}*[-{a1q1*a1} - {a2q2*a2} + {a2q2*a1} + {a1q1*a2}] = {_dE*self.beta}')
                print(f'              _p = {min(1.0, np.exp(-_dE*self.beta))}')
        else:
            e1 = struc1.getEnergy()
            e2 = struc2.getEnergy()
            dE = -e1*a1 -e2*a2 + e2*a1 + e1*a2 # delta = (e1-e2)*(alp2-alp1)
        if np.isnan(dE):
            return float('NaN')
        return min(1.0, np.exp(-dE*self.beta))

    def _exchangeXyz(self, w1, w2):
        _, struc1, sim1, *_ = w1
        _, struc2, sim2, *_ = w2
        xyz1 = struc1.getXyz()
        xyz2 = struc2.getXyz()
        e1 = struc1.getEnergy()
        e2 = struc2.getEnergy()
        struc1.changeXyz(xyz2, e2)
        struc2.changeXyz(xyz1, e1)
        rf1 = sim1.getRestartFile()
        rf2 = sim2.getRestartFile()
        sim1.updateRestartFile(rf2)
        sim2.updateRestartFile(rf1)
        sim1.updateXyz(xyz2)
        sim2.updateXyz(xyz1)

    def _MHCycle(self):
        skipNext = False
        order = range(len(self.alphaSet)-1)
        if self.changeOrder and self.Niter%2!=0:
            order = range(len(self.alphaSet)-2,-1,-1)
        for i in order:
            if skipNext:
                skipNext = False
                continue
            w1 = self.alphaSet[i]
            w2 = self.alphaSet[i+1]
            a1, _, sim1, *_ = w1
            a2, _, sim2, *_ = w2
            if (sim1.results['exitCode']|sim2.results['exitCode']):
                print(f'skip exchange attempt of {a1} {a2}: non-zero exit code')
                continue
            p0 = self.randomGenerator.random()
            p = self._calcP(w1,w2)
            # p=1 # REMOVE
            if p>=p0:
                indices = (i,i+1)
                pair = (w1,w2,p)
                if self.unrestrictedExchange:
                    for j in range(i+2,len(self.alphaSet)):
                        w3 = self.alphaSet[j]
                        a3, _, sim3, *_ = w3
                        if (sim3.results['exitCode']):
                            print(f'skip exchange attempt of {a1} {a3}: non-zero exit code')
                            continue
                        p0 = self.randomGenerator.random()
                        p = self._calcP(w1,w3)
                        if p>=p0:
                            indices = (i,j)
                            pair = (w1,w3,p)
                        else:
                            break
                self._exchangeXyz(*pair[:-1])
                self.Counter.updateIndex(*indices)
                # print(f'Pair {pair[0][0]}[{indices[0]}] and {pair[1][0]}[{indices[1]}] exchanged with p = {pair[2]}')
                print(f'Pair {pair[0][0]}[{indices[0]}] and {pair[1][0]}[{indices[1]}] exchanged with p = {pair[2]:.3f}')
                self.swapCount[i]+=1
                skipNext = True

    def runREMD(self, steps):
        for step in range(1, steps):
            print(f'Step {step}. Start MD calculations\n')
            self._runMD()
            print(f'\tTry to exchange')
            self._MHCycle()
            if step%self.Ncheck==0:
                if self.addWorlds and len(self.alphaSet) < self.Nmax:
                    alphas = [w[0] for w in self.alphaSet]
                    offset = 0
                    for i,n in enumerate(self.swapCount):
                        if n <= self.Ncheck*self.Pmin//100: # т.е., обмены данной пары происходили с вероятностью менее Pmin
                            a1 = alphas[i]
                            a2 = alphas[i+1]
                            self.addWorld(0.5*(a1+a2), i+1+offset, parm = self.initParm, vars = self.initVars, options = self.initOptions, *self.initOther)
                            offset+=1
                if self.delWorlds:
                    for i,n in enumerate(self.swapCount[1:]):
                        if n >= self.Ncheck*self.Pmax//100:
                            self.removeWorld(i)
                self.swapCount = [0 for _ in range(len(self.alphaSet)-1)]
            self.Counter.addRow()
            self.writeStat()
            self.Niter+=1
            print(f'Iteration {self.Niter} done')
        self.Counter.writeLog()
        print(f'\nresulting set of worlds:\n\n{[w[0] for w in self.alphaSet]}\n')
        print('aREMD calculation done')

    def addWorld(self, alpha, i, wd = None, datafile = None, parm = None, vars = None, options = None, *other):
        assert 0.0<=alpha<=1.0, f'invalid alpha value: {alpha}'
        if wd == None:
            wd = os.path.join(self.baseDir, f'{self.Nadd}-{alpha:.3e}')
        if datafile == None:
            datafile = self.initDataFile
        sim = engines.Simulation(alpha, wd, datafile, parm=parm, vars=vars, options=options)
        ps = self.initStruc
        xyz = ps.getXyz()
        e = ps.getEnergy()
        struc = Structure(ps.topology, ps.a, ps.b, ps.c, ps.alpha, ps.beta, ps.gamma, xyz=xyz) # TODO: make it more clear
        struc.e = e
        self.alphaSet.insert(i,(alpha, struc, sim, *other))
        self.Counter.addIndex(i)
        self.Nadd+=1
        print(f'New world with alpha = {alpha} was added')

    def removeWorld(self, i):
        assert 0 < i < len(self.alphaSet)-1, 'invalid index'
        alpha = self.alphaSet[i][0]
        self.alphaSet.pop(i)
        self.Counter.delIndex(i)
        self.Ndel+=1
        print(f'World with alpha = {alpha} was deleted')

    def writeStat(self):
        energy = []
        for world in self.alphaSet:
            _, struc, *_ = world
            energy.append(struc.getEnergy())
        with open(self.dumpEnergy,'a') as fo:
            fo.write(''.join([f'{e:.6e}'.ljust(20) for e in energy]) + '\n')
        # xyz = struc.getXyz()
        # e = energy[-1]
        # lines = xyz.split('\n')
        for i,world in enumerate(self.alphaSet):
            if i in self.idump:
                alpha, struc, *_ = world
                e = struc.getEnergy()
                xyz = struc.getXyz()
                print(f'dump xyz from H({alpha:.6f})')
                xyz[1]=f'Energy = {e:.2f} kcal/mol\n'
                with open(self.dumpXYZ,'a') as fo:
                    for row in xyz:
                        fo.write(row)
                if i==len(self.alphaSet)-1:
                    with open(self.dumpXYZ[:-4]+'_w1.xyz','a') as fo:
                        for row in xyz:
                            fo.write(row)

def optimize(args, kwargs):
    wd, datafile, xyz, cores,indx = args
    try:
        os.mkdir(wd)
    except:
        shutil.rmtree(wd)
        os.mkdir(wd)
    # print(kwargs)
    # print(*kwargs)
    sim = engines.Simulation(1.0, wd, datafile, xyz = xyz, **kwargs)
    args = (wd, sim.IN_FNAME, sim.vars['d'], sim.ERR_FNAME, sim.TRJ_FNAME, sim.k, sim.patt, cores)
    sim.prepare()
    e0, e1, xyzopt = engines.runOpt(args)
    print(f'[index {indx}] Optimization: {e0:.2f} --> {e1:.2f} kcal/mol')
    if xyzopt != None:
        xyzopt[1] = f'Energy = {e1:.2f} kcal/mol\n'
    shutil.rmtree(wd)
    return xyzopt,indx

def optimizeFrames(trjfile, datafile, sort=True, maxp = 12, cores=1, **kwargs):
    args = []
    XYZframes = base_utils.readXYZ(trjfile)
    wd = os.path.dirname(os.path.abspath(trjfile))
    optdir = os.path.join(wd, f'opt-{time.time()}')
    os.mkdir(optdir)
    os.chdir(optdir)
    for i,frame in enumerate(XYZframes):
        fname = f'{i+1}.xyz'
        with open(fname,'w') as fo:
            for line in frame:
                fo.write(line)
        args.append([(os.path.join(optdir, f'opt-{i+1}'), datafile, os.path.abspath(fname), cores,i),kwargs])
    # print(args)
    with mp.Pool(maxp) as pool:
        # opt_frames = pool.starmap(optimize, args)
        opt_frames = pool.starmap(optimize, args)
    opt_frames = [i for i in opt_frames if i[0]!=None]
    best_frame,indx = min(opt_frames, key= lambda x: float(x[0][1].split()[2]))
    print(f'best structure with E = {float(best_frame[1].split()[2])} with index = {indx}')
    with open(os.path.join(wd,trjfile[:-4]+'-best.xyz'),'w') as fo:
        for line in best_frame:
            fo.write(line)
    
    if sort:
        opt_frames = sorted(opt_frames, key= lambda x: float(x[0][1].split()[2]))
    os.chdir(wd)
    optfxyz = os.path.join(wd,trjfile[:-4]+'-opt.xyz')
    with open(optfxyz,'w') as fo:
        for frame in opt_frames:
            xyz,indx = frame
            for line in xyz:
                fo.write(line)
    print(f'Optimization of {len(XYZframes)} frames done')
    print(f'Number of successful jobs: {len(opt_frames)}')
    print(f'Number of failed jobs: {len(XYZframes)-len(opt_frames)}')
    shutil.rmtree(optdir)
    # write_data(wd,datafile,optfxyz,kwargs)
    # print('data file with relaxed structure: %s' % optfxyz)



def write_data(wd, data, fxyz,**kwargs):
    # with tempfile.NamedTemporaryFile(dir=wd) as fi:
    sim = engines.Simulation(1.0, wd, data,**kwargs)
    _xyz = base_utils.readXYZ(fxyz, 0)
    sim.updateXyz(_xyz)
    sim.prepare(options={'blank':True})
    sim.DAT_FNAME = 'REMD_opt.data'
    sim.TRJ_FNAME = ''
    args = (sim.WD, sim.IN_FNAME, sim.vars['d'], sim.ERR_FNAME, sim.TRJ_FNAME, sim.k, sim.patt, 1)
    res = engines.runOpt(args)
    sim.update(res)



# test for single point
# elem = base_utils.getElementsLmp('/home/md/md/HREMD_v1/toy.lmps')
# struc1 = Structure('/home/md/md/HREMD_v1/toy.lmps','','','','','','','')
# struc2 = Structure('/home/md/md/HREMD_v1/toy.lmps','','','','','','','')
# sim1 = engines.Simulation(1.0,'/home/md/md/HREMD_v1/single_point_test','/home/md/md/HREMD_v1/toy.lmps', parm = {'improper_style' : 'umbrella', 'dihedral_style' : 'harmonic','dump_modify' : 'DUMPFILE element '+' '.join(elem)})
# sim2 = engines.Simulation(0.01,'/home/md/md/HREMD_v1/single_point_test','/home/md/md/HREMD_v1/toy.lmps', parm = {'improper_style' : 'umbrella', 'dihedral_style' : 'harmonic','dump_modify' : 'DUMPFILE element '+' '.join(elem)})
# world1 = (1.0, struc1, sim1)
# world2 = (0.01, struc2, sim2)
# w1s1=REMD.calc_single_point(world1,struc1)
# w2s2=REMD.calc_single_point(world2,struc2)
# w1s2=REMD.calc_single_point(world1,struc2)
# w2s1=REMD.calc_single_point(world2,struc1)
# print(f'dE = {-w1s1-w2s2+w2s1+w1s2}')
# exit()



# global parms: T, N, seed, units, LAMMPS path
# global options: [options] + doMTD, unrestrictedExchange, selectLastStruc, addWorlds, delWorlds + Pmin, Pmax, Ncheck
if __name__=='__main__':
    parm = {'improper_style' : 'umbrella', 'dihedral_style' : 'harmonic'}
    write_data('/home/md/md/PES-185_dens/08.09-TEST','/home/md/pysimm/pysimm/PES_185_dens_box/x8_90_mono/uniform_polymer_soft.lmps','/home/md/md/PES-185_dens/08.09-TEST/structures-opt.xyz', parm=parm)
    exit()
    # alphaRange = [0.7, 0.8, 0.9, 1.0]
    alphaRange = [0.0, 1.0]
    # initStruc, alphaRange, datafile, baseDir, parm=None, vars=None, options=None, *other

    # topology,xyz, a, b, c, alpha, beta, gamma, symm = 'P1'
    # xyz1 = base_utils.get_frame_xyz('/home/artem/LAMMPS_TEST/traj.xyz')
    # xyz2 = base_utils.getXyzfromData('/home/artem/HREMD/src/lammps.data')
    # spl1 = xyz1.split('\n')
    # spl2 = xyz2.split('\n')
    # print(len(spl1))
    # print(len(spl2))
    # with open('fromTraj.xyz','w') as f1:
    #     f1.write(''.join(xyz1))
    # with open('fromData.xyz','w') as f2:
    #     f2.write(''.join(xyz2))


    dataFile = '/home/artem/HREMD/src/lammps.data'
    istruc = Structure(dataFile,None,None,None,None,None,None)
    istruc.e = float('nan')#-1.807965e+02

    baseDir = '/home/artem/LAMMPS_TEST/dbg-10.04/'
    options = {
        'minBeforeMD' : True,
        'doLongES' : True
    }
    Set1 = REMD(istruc, alphaRange, dataFile, baseDir)

    Set1.runREMD(40)

    print([w[0] for w in Set1.alphaSet])
    print([str(w[0]) for w in Set1.alphaSet])

