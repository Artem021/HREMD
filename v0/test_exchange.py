import random, math

random.seed(228)

test_add = True
test_del = False


EXCHANGE_STAT = 'excDBG.data-REF'

# alpha_set = [0.0, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
alpha_set = [0.0, 1.0]

open(EXCHANGE_STAT,'w').close()

indices = [i for i in range(len(alpha_set))]
imap = {i:indices.index(i) for i in indices}

rows = []
deleted = []

def addIndex(indices, rows, j):
    # добавление происходит справа: т.е. i-й элемент < j остается при своем значении, 
    # j становится j, а элементы после j смещаются на 1; при этом j не должен замещать собой первый и последний элементы 
    assert 0 < j <= len(indices)-1
    for i in range(len(indices)):
        if indices[i]>=j:
            indices[i]+=1
    for row in rows:
        for i in range(len(indices)):
            if row[i]>=j:
                row[i]+=1
        row.insert(j,float('Nan'))
    indices.insert(j,j)
    imap = {i:i for i in range(len(indices))}

    return indices, imap, rows

def delIndex(indices, deleted, j):
    assert 0 < j < len(indices)-1
    val = indices.pop(j)
    for i in range(len(indices)):
        if indices[i]>val:
            indices[i]-=1
    imap = {i:i for i in range(len(indices))}
    deleted.append(j)
    return indices, imap, deleted

def updateIndex(indices, i):
    i1 = indices.index(i)
    i2 = indices.index(i+1)
    atmp = indices[i1]
    indices[i1] = indices[i2]
    indices[i2] = atmp
    return indices

def addRow(imap, indices, deleted, rows):
    row = [imap[i] for i in indices]
    for i in reversed(deleted):
        for j in range(len(row)):
            if row[j]>=i:
                row[j]+=1
        row.insert(i, float('Nan'))
    rows.append(row)
    return row, rows

Nmax=400
niter=1
exc = False
while niter<Nmax:
    print('-'*50+f'\nIteration {niter}\n'+'-'*50)
    for i in range(len(alpha_set)-1):
        alp1 = alpha_set[i]
        alp2 = alpha_set[i+1]
        if exc:
            exc=False
            print(f'Skip : {alp1} : {alp2}')
            continue
        print(f'Attempt : {alp1} : {alp2}')
        p1 = random.random()
        p2 = math.exp(-50*abs((alp1-alp2)))
        if p2>p1:
            indices = updateIndex(indices, i)
            print(f'exchange: {alp1} --> {alp2}')
            exc=True
    if niter%20==0:
        if test_add:
            newI = random.choice([i for i in range(len(alpha_set))][1:])
            newA = 0.5*(alpha_set[newI] + alpha_set[newI-1])
            alpha_set.insert(newI,newA)
            indices, imap, rows = addIndex(indices, rows, newI)
            print(f'new world with a={newA} was created')
        if test_del:
            p = random.random()
            if p<=0.33:
                newI = random.choice([i for i in range(len(alpha_set))][1:-1])
                oldA = alpha_set.pop(newI)
                indices, imap, deleted = delIndex(indices, deleted, newI)
                print(f'world with a={oldA} was deleted')

    row, rows = addRow(imap, indices, deleted, rows)
    niter+=1
    if len(row)==15:
        test_add=False
        test_del=False



print(alpha_set)
print(len(rows))
print(niter)
with open(EXCHANGE_STAT,'a') as exc_stat:
    for row in rows:
        exc_stat.write('   '.join([str(i) for i in row]) + '\n')
        # exc_stat.write('   '.join([str(exchange_ord[a]) for a in exchange_set]) + '\n')
        # exc_stat.write('   '.join([str(alpha_ord[i]) for i in alpha_indices]) + '\n')
        # exc_stat.write('   '.join([str(alpha_indices[i]) for i in sorted(alpha_indices)]) + '\n')
