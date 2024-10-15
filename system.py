from pysimm import system, lmps, forcefield
from pysimm.apps.random_walk import random_walk, copolymer
from multiprocessing import Pool
CPU_PER_TASK = 24
MAX_TASK = 10


# monomer A (Bisphenol A + 4,4'-Dichlorodiphenyl sulfone)
monomerA = '/home/md/pysimm/pysimm/PES_185_50/topo/PES_185.mol'
ih = 11
it = 50


def monomer(molfile,ih,it,position):
    sys = system.read_mol(molfile)
    ff = forcefield.Dreiding()
    sys.apply_forcefield(ff)
    head = sys.particles[ih]
    tail = sys.particles[it]
    head.linker = 'head'
    tail.linker = 'tail'
    if position=='tail' or position=='middle':
        for bond in head.bonds: # tail
            if bond.a.elem == 'H' or bond.b.elem == 'H':
                hatom = bond.a if bond.b is tail else bond.b
                sys.particles.remove(hatom.tag, update=False)
                break
    if position=='head' or position=='middle':  
        for bond in tail.bonds: # head 
            if bond.a.elem == 'H' or bond.b.elem == 'H':
                hatom = bond.a if bond.b is tail else bond.b
                sys.particles.remove(hatom.tag, update=False)
                break
    sys.remove_spare_bonding()
    sys.pair_style = 'lj/cut'
    lmps.quick_min(sys, min_style='fire')
    sys.add_particle_bonding()
    return sys


from pysimm import system, lmps, forcefield
from pysimm.apps.random_walk import random_walk
# from midle import monomer as monomerM
# from tail import monomer as monomerT
# from head import monomer as monomerH
from pysimm.apps.random_walk import copolymer



# TODO: expand to copolymers

def build_chain(size, monomers : list[system.System], ratio : list[float] = None):
    chainM = monomer()
    chainH = monomer()
    chainT = monomer()
    chainM.pair_style = 'lj'
    chainH.pair_style = 'lj'
    chainT.pair_style = 'lj'
    f = forcefield.Dreiding()
    chainM.apply_charges(f, charges='gasteiger')
    chainH.apply_charges(f, charges='gasteiger')
    chainT.apply_charges(f, charges='gasteiger')
    # chain = copolymer([chainT, chainM, chainH], size, pattern=[1, size-2, 1], forcefield=f,settings={'np':CPU_PER_TASK,'prefix':'mpirun'})
    chain = copolymer([chainT, chainM, chainH], size, pattern=[1, size-2, 1], forcefield=f)
    chain.apply_forcefield(f)
    chain.apply_charges(f, charges='gasteiger')
    return chain


"""def build_chain(size, monomers):
    chainM = monomer()
    chainH = monomer()
    chainT = monomer()
    chainM.pair_style = 'lj'
    chainH.pair_style = 'lj'
    chainT.pair_style = 'lj'
    f = forcefield.Dreiding()
    chainM.apply_charges(f, charges='gasteiger')
    chainH.apply_charges(f, charges='gasteiger')
    chainT.apply_charges(f, charges='gasteiger')
    # chain = copolymer([chainT, chainM, chainH], size, pattern=[1, size-2, 1], forcefield=f,settings={'np':CPU_PER_TASK,'prefix':'mpirun'})
    chain = copolymer([chainT, chainM, chainH], size, pattern=[1, size-2, 1], forcefield=f)
    chain.apply_forcefield(f)
    chain.apply_charges(f, charges='gasteiger')
    return chain"""


def build_chains_mp(sizes):
    with Pool(MAX_TASK) as pool:
        return pool.map(build_chain, sizes)
        
def build_chains(sizes):
    return [build_chain(n) for n in sizes]

def build(chains, dens=1.1,name='uniform_polymer'):
    uniform_polymer = system.replicate(chains, [1 for i in chains], density=dens, rand=True)
    uniform_polymer.write_xyz(f'{name}.xyz')
    uniform_polymer.write_yaml(f'{name}.yaml')
    uniform_polymer.write_lammps(f'{name}.lmps')
    uniform_polymer.write_chemdoodle_json(f'{name}.json')



def run(test=False):
    # we'll create a pe monomer from the pysimm.models database
    chainM = monomerM()
    chainH = monomerH()
    chainT = monomerT()

    chainM.pair_style = 'lj'
    chainH.pair_style = 'lj'
    chainT.pair_style = 'lj'
    
    f = forcefield.Dreiding()
    chainM.apply_charges(f, charges='gasteiger')
    chainH.apply_charges(f, charges='gasteiger')
    chainT.apply_charges(f, charges='gasteiger')      
    
    # chainT.write_xyz('chainT.xyz')
    # chainT.write_lammps('chainT.lmps')
    # chainM.write_xyz('chainM.xyz')
    # chainM.write_lammps('chainM.lmps')
    # chainH.write_xyz('chainH.xyz')
    # chainH.write_lammps('chainH.lmps')
    # exit()
    
    

    print('Building polymer chain 1...')
    # monochain = copolymer([chainT, chainM, chainH], 82, pattern=[1, 80, 1], forcefield=f)
    monochain = copolymer([chainT, chainM, chainH], 20, pattern=[1, 18, 1], forcefield=f) # 2.16, 82
    # monochain = copolymer([chainT, chainM, chainH], 14, pattern=[1, 12, 1], forcefield=f) # 2.16, 82
    monochain.apply_forcefield(f)
    monochain.apply_charges(f, charges='gasteiger')
    monochain.write_xyz('monochain.xyz')
    monochain.write_yaml('monochain.yaml')
    monochain.write_lammps('monochain.lmps')
    monochain.write_chemdoodle_json('monochain.json')
    
    # print('Building polymer chain 2...')
    # # monochain2 = copolymer([chainT, chainM, chainH], 82, pattern=[1, 80, 1], forcefield=f)
    # monochain2 = copolymer([chainT, chainM, chainH], 4, pattern=[1, 2, 1], forcefield=f) # 2.16, 82
    # monochain2.apply_forcefield(f)
    # monochain2.apply_charges(f, charges='gasteiger')
    # monochain2.write_xyz('monochain2.xyz')
    # monochain2.write_yaml('monochain2.yaml')
    # monochain2.write_lammps('monochain2.lmps')
    # monochain2.write_chemdoodle_json('monochain2.json')
    
    
    
    
    
    
    
    print('Replicating polymer chain...')
    uniform_polymer = system.replicate(monochain, 8, density=1.1, rand=True)
    # uniform_polymer = system.replicate([monochain,monochain2], [1,1], density=1.1, rand=True)
    uniform_polymer.write_xyz('uniform_polymer.xyz')
    uniform_polymer.write_yaml('uniform_polymer.yaml')
    uniform_polymer.write_lammps('uniform_polymer.lmps')
    uniform_polymer.write_chemdoodle_json('uniform_polymer.json')
    
if __name__ == '__main__':
    # run()
    import time
    t0=time.time()
    chains = build_chains([46, 71, 49, 51, 56, 138, 21, 48])
    build(chains,dens=1.1,name='uniform_polymer')
    t1=time.time()
    print(f'done in {int(t1-t0)} sec')
    