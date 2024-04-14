import numpy as np, os

LOG = '/home/artem/HREMD_algorithm/cp2k/slurm-115895.out'


pairs = {}
worlds = {}

read_set = []
with open(LOG,'r') as fi:
    while True:
        line = fi.readline()
        if line=='':
            break
        if 'exchange attempt' in line:
            if '0.0 and 0.001' in line:
                read_set = []
            w1,_,w2 = line.split()[-3:]
            pair = '_'.join([w1,_,w2])

            iline = fi.readline()
            assert 'first index' in iline
            iline = iline.replace(';','')
            i1,_,_,i2 = iline.split()[-4:]

            if w1 not in read_set:
                try:
                    worlds[w1].append(i1)
                except:
                    worlds[w1] = [i1]
                read_set.append(w1)
            if w2 not in read_set:
                try:
                    worlds[w2].append(i2)
                except:
                    worlds[w2] = [i2]
                read_set.append(w2)

            eline = fi.readline()
            assert 'Relative energy and Delta' in eline
            rE,_,delta = eline.split()[-3:]
            
            pline = fi.readline()
            assert 'p and p0' in pline
            p,_,p0 = pline.split()[-3:]

            excline = fi.readline()
            assert 'exchange:' in excline
            exc = excline.split()[-1]
            exc = 'True' in exc

            try:
                pairs[pair]['rE'].append(rE)
                pairs[pair]['delta'].append(delta)
                pairs[pair]['p'].append(p)
                pairs[pair]['p0'].append(p0)
                pairs[pair]['exc'].append(exc)
            except:
                pairs[pair] = {
                    'rE' : [],
                    'delta': [],
                    'p' : [],
                    'p0' : [],
                    'exc' : []
                }

os.chdir('/home/artem/HREMD_algorithm/PAM-crystal/PAM-12ns/')
skeys_w = sorted(worlds.keys(), key=lambda x: float(x))
skeys_p = sorted(pairs.keys(), key=lambda x: min([float(a) for a in x.split('_and_')]))

with open('lowest_frames.dat','w') as fo:
    n = 0
    nmin = min(len(worlds[x]) for x in worlds)
    nmax = max(len(worlds[x]) for x in worlds)
    for x in worlds:
        print(x, len(worlds[x]))
    assert nmin==nmax
    # head = '               '.join(skeys_w) + '\n'
    head = [f'{x}'.center(15) for x in skeys_w]
    fo.write(''.join(head) + '\n')
    while n < nmin:
        raw = [f'{worlds[x][n]}'.center(15) for x in skeys_w]
        line = ''.join(raw) + '\n'
        fo.write(line)
        n+=1

ncenter = 30
with open('exchange_stat.dat','w') as fo:
    n = 0
    nmin = min(len(pairs[pair][val]) for val in pairs[pair] for pair in pairs)
    nmax = max(len(pairs[pair][val]) for val in pairs[pair] for pair in pairs)
    assert nmin==nmax
    head = [f'{x}'.center(ncenter) for x in skeys_p]
    fo.write('Property'.center(ncenter)+''.join(head) + '\n')
    while n < nmin:
        for val in pairs[skeys_p[0]]:
            raw = [f'{float(pairs[p][val][n]):.6f}'.center(ncenter) for p in skeys_p]
            line = val.center(ncenter) + ''.join(raw) + '\n'
            fo.write(line)
        n+=1
print(0)

