import os


WD = '/home/artem/PAM/24.10-6ns/'

diagram = 'exchanges.log'

worldX = '0'
I = 0

alpha_set = {}

dir_set = []


def read_exchanges(log):
    M = []
    n = 0
    with open(log,'r') as fi:
        for line in fi:
            worlds = [int(i) for i in line.split()]
            M.append(worlds)
            n+=1
    return n, M



def merge_trajectory():
    pass


os.chdir(WD)

for d in os.listdir():
    if os.path.isdir(d):
        try:
            print(float(d))
            dir_set.append(d)
        except:
            continue

dir_set = sorted(dir_set, key=lambda a: float(a))
nmax, M = read_exchanges(diagram)
assert len(M[0])==len(dir_set)
alpha_set = {i:dir_set[i] for i in range(len(dir_set))}


n=1
with open(f'world-{worldX}-evolution.xyz','w') as traj:
    for row in M[:10]:
        index = row[I]
        with open(f'{alpha_set[index]}/PAM-pos-1-{n}.xyz') as fi:
            traj.write(fi.read())
            print(f'add fragment from {alpha_set[index]}/PAM-pos-1-{n}.xyz')


print(alpha_set)
print(dir_set)


