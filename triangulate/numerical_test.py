import math
import numpy as np
import scipy
import scipy.optimize
import csv
from scipy.optimize import fsolve
import sympy as sp

def get_vol_idx(heights):
    h1 = heights[0]
    h2 = heights[1]
    h3 = heights[2]
    idx_1 = 0
    idx_2 = 0
    idx_3 = 0
    for i in range(1,h1+1):
        idx_1 += i*(i+1)/2
    idx_1 -= 1
    idx_2 = (h2)*(h2+1)/2
    idx_3 = h3
    idx = idx_1 + idx_2 + idx_3
    return idx

def func(p, *d):
    f1, f2, f3 = d
    return(eval(f1), eval(f2), eval(f3))
    #return (f1.evalf(b1 = p[0], b2 = p[1], b3 = p[2]), f2.evalf(b1 = p[0], b2 = p[1], b3 = p[2]), f3.evalf(b1 = p[0], b2 = p[1], b3 = p[2]))
    #return (f1.evalf(subs={b1 : p[0], b2 : p[1], b3 : p[2]}), f2.evalf(subs={b1 : p[0], b2 : p[1], b3 : p[2]}), f3.evalf(subs={b1 : p[0], b2 : p[1], b3 : p[2]}))
    
def constraint(Series, sol):
    global vol_min_global
    vol = Series.evalf(subs={b1 : sol[0], b2 : sol[1], b3 : sol[2]})
    vol = abs(vol)
    if vol <= 1 and vol >= vol_min_global:
        return 1, vol
    
    print ('volume: ', vol, ' is out of bounds.')
    
    return 0, -1

def test_func_1(x, a):
    return a * float(1.0)/(x)
def test_func_2(x, a, b, c):
    return -1*a * np.log(x-b) + c
def test_func(x):
    return test_func_1((x+10), 4.62386348)+0.03
def test_func_inv(y):
    return 4.62386348/(y-0.03)-10
def dist_to_func(pt):
    x0 = float(pt[0])
    y0 = float(pt[1])
    return abs(y0 - test_func(x0))

def NSolve(Series, d, bound):
    vol = -1
    sol = -1
    const = 0
    count = 0
    MAX_COUNT = 3
    
    b1_min = bound[0][0]
    b1_max = bound[0][1]
    b2_min = bound[1][0]
    b2_max = bound[1][1]
    b3_min = bound[2][0]
    b3_max = bound[2][1]

    while const == 0:
        if count >= MAX_COUNT:
            return vol,sol
            
        count += 1
        d1_0 = np.random.uniform(low=b1_min, high=b1_max)
        d2_0 = np.random.uniform(low=b2_min, high=b2_max)
        d3_0 = np.random.uniform(low=b3_min, high=b3_max)
        print ('reset starting point: ', d1_0, d2_0, d3_0)

#         try:
        sol = fsolve(func, x0 = np.array([d1_0, d2_0, d3_0]), args = d)
        print ('solution: ', sol)
        print ('guessed vol: ', Series.evalf(subs={b1 : sol[0], b2 : sol[1], b3 : sol[2]}))
#         except:
#             continue
        
        const, vol = constraint(Series, sol)

    print ('Done.')

    return vol, sol

def fit_NSolve(Series, max_range_start, heights, vol_range):
    d1 = sp.diff(Series, b1)
    d2 = sp.diff(Series, b2)
    d3 = sp.diff(Series, b3)
    d1 = str(d1).replace('b1', 'p[0]').replace('b2', 'p[1]').replace('b3', 'p[2]')
    d2 = str(d2).replace('b1', 'p[0]').replace('b2', 'p[1]').replace('b3', 'p[2]')
    d3 = str(d3).replace('b1', 'p[0]').replace('b2', 'p[1]').replace('b3', 'p[2]')
    d = (d1, d2, d3)
    MAX_NUM_VOL = 10
    target_vol = vol_range[0]
    max_diff = vol_range[1]
    range_dict = {}
    height = max(heights)
    idx = get_vol_idx(heights)
    print ('vol idx: ', idx)
    
    fit_dist_min = 999
    vol_ret = -1
    sol_ret = -1
    
    for max_range in range(max_range_start, max_range_start + height):
        vol_list = []
        sol_list = []
        dist_list = []
        for i in range(1, max_range+1):
            if len(vol_list) > MAX_NUM_VOL:
                break
            for j in range(1, max_range+1):
                for k in range(1, max_range+1):
                    try:
                        bound = range_dict['%d_%d_%d' % (i,j,k)]
                        continue
                    except:
                        pass
                    bound = [[i-1,i],[j-1,j],[k-1,k]]
                    range_dict['%d_%d_%d' % (i,j,k)] = bound
                    print ('try bound: ', bound)
                    range_dict['%d_%d_%d' % (i,j,k)] = bound
                    vol, sol = NSolve(Series, d, bound)
                    if type(sol) == int or type(vol) == int:
                        print ('range ', i,j,k,' does not work')
                        continue
                    if type(sol) == np.ndarray:
                        sol = sol.tolist()
                    target_dist = target_vol - vol
                    fit_dist = dist_to_func([idx, vol])
                    print ('dist to fit function: ', fit_dist)
                    if abs(fit_dist) < 1e-6:
                        return vol, sol
                    if idx < 100:
                        if vol <= target_dist:
                            return vol, sol
                        else:
                            vol_list.append(vol)
                            sol_list.append(sol)
                            dist_list.append(fit_dist)
                    else:
                        if fit_dist < max_diff:
                            if vol <= target_dist:
                                return vol, sol
                            else:
                                vol_list.append(vol)
                                sol_list.append(sol)
                                dist_list.append(fit_dist)
        print ('Done.')
        vol_list = list(vol_list)
        sol_list = list(sol_list)
        print ('vol_list: ', vol_list)
        print ('sol_list: ', sol_list)
        print ('dist_list: ', dist_list)
        if len(vol_list) != 0 and len(sol_list) != 0:
            for i in range(len(vol_list)):
                vol_tmp = vol_list[i]
                sol_tmp = sol_list[i]
                dist = dist_list[i]
                if abs(dist) < 1e-6:
                    return vol_tmp, sol_tmp
                if abs(dist) < max_diff:
                    if abs(dist) < fit_dist_min:
                        fit_dist_min = abs(dist)
                        vol_ret = vol_tmp
                        sol_ret = sol_tmp
                    continue
                else:
                    print ('distance ', dist, ' is too far from target')
                    continue
        try:
            if type(sol_ret) != int and len(sol_ret) != 0:
                print ('good')
                return vol_ret, sol_ret
        except:
            continue
    print ('no valid solution')
    return vol_ret, sol_ret

def parse_series(height_file, input_path_list):
    b1 = sp.symbols('b1')
    b2 = sp.symbols('b2')
    b3 = sp.symbols('b3')
    start = 1j
    end = 1j
    path_idx = 0
    for line in height_file:
        h = eval(line)
        h_start = h[0]
        h_end = h[1]
    print (h_start, h_end)
    for input_path in input_path_list:
        counter = 0
        try:
            input_file = open(input_path, 'r')
        except:
            continue
        for line in input_file:
            data_idx = line.find(']', 2)
            h_list = list(eval(line[2:data_idx]))
            #print (h_list, ' idx: ', counter)
            if h_list == h_start:
                start = [counter, path_idx]
                print ('start: ', start)
            elif h_list == h_end:
                end = [counter, path_idx]
                print ('end: ', end)
                input_file.close()
                return start, end
            counter += 1
        path_idx += 1
        input_file.close()
    print ('end for real')
    return start, end

def find_vol(input_path_list, start, end, out_path):
    global vol_min_global
    toggle = 0
    path_idx = 0
    for input_path in input_path_list:
        counter = 0
        try:
            input_file = open(input_path, 'r')
        except:
            continue
        for line in input_file:
            data_idx = line.find(']', 2)
            h_list = list(eval(line[2:data_idx]))
            if counter == start[0]:
                print ('Start')
                toggle = 1
            if toggle == 1:
                print ('counter: ', counter)
                data = eval(line)
                h = data[0]
                series = data[1]
                h_max = max(h)
                print (h, 16.0/27/h_max)
                vol_min_global = 1.0/(h_max**3)
                vol, sol = fit_NSolve(series, 3, h, [16.0/27/h_max, 0.02])
                print ('vol: ', vol, ' sol: ', sol)
                out_file = open(out_path, 'a')
                out_file.write('[%s, %f, %s]\n'%(str(data[0]), vol, str(sol)))
                out_file.close()
                print (vol, sol)
                if counter == end[0] and path_idx == end[1]:
                    input_file.close()
                    print ('Done.')
                    return
            counter += 1
        path_idx += 1
        input_file.close()
        print ('Finished.')

h_min = 32
h_max = 32
h_idx = 1
input_path_0 = 'series_%d_%d.txt'%(h_min, h_max)
input_path_1 = 'series_%d_%d.txt'%(h_min+1, h_max+1)
input_path_2 = 'series_%d_%d.txt'%(h_min+2, h_max+2)
input_path_list = [input_path_0, input_path_1, input_path_2]
out_path = 'test/vol_%d_%d.txt'%(h_min, h_max)
out_file = open(out_path, 'w')
out_file.close()
vol_min_global = 1/2.0**3
b1 = sp.symbols('b1')
b2 = sp.symbols('b2')
b3 = sp.symbols('b3')
height_path = 'series/height_%d.txt' % h_idx
height_file = open(height_path, 'r')
start, end = parse_series(height_file, input_path_list)
print (start, end)
find_vol(input_path_list, start, end, out_path)