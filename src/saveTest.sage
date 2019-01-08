from pathlib import Path

import numpy as np
import os

def write_file(fname, a):
    p = Path(fname)
    test_file = p.open('ab')
    np.save(test_file, a)
    test_file.close()

def load_file(fname):
    p = Path(fname)
    out_list = []
    with p.open('rb') as f:
        fsz = os.fstat(f.fileno()).st_size
        out = np.load(f)
        out_list.append(out)
        while f.tell() < fsz:
            out_list.append(np.load(f))
    return out_list

triang = load_file('triang.npy')
power = load_file('power.npy')

print('Number of triang: ',len(triang)) #6558
print('Number of power: ', len(power))

data_file = open('lift_1_40.txt', 'r')
data_x = []
for line in data_file:
    data = eval(line)
    if data[1] > 0:
        data_x.append(data[0])
for i in range(100):      
    print(i,'-th heights:')
    print(data_x[i])
    
    print(i,'-th triang:')
    print(triang[i][-1])

print('Done.')
