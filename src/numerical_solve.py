#/usr/bin/python
"""
:mod:`numerical solver` -- Find the minimum volume of a polytope
===================================
 
.. module:: numerical_solver
   :platform: Linux (Ubuntu 16.04)
   :synopsis: highlight document snippets that match a query.
.. moduleauthor:: Antonio Cao (antonio.cao@yale.edu)
 
  
Requirements::
    1.  You will need to install the numpy and scipy library to run this code.
        https://scipy.org/install.html
 
"""
import numpy as np
import sympy as sp
from sympy import *
from sympy.utilities.lambdify import lambdify
import scipy.optimize
import csv
from scipy.optimize import fsolve
import os
from os.path import expanduser
from pathlib import Path
import itertools

def merge_lists(list1, list2):
    return [x for x in itertools.chain.from_iterable(itertools.izip_longest(list1,list2)) if x]

def find_bound(i,j,k, bound):
    b1_min = bound[0][i]
    b1_max = b1_min+1
    b2_min = bound[1][j]
    b2_max = b2_min+1
    b3_min = bound[2][k]
    b3_max = b3_min+1
    return [[b1_min, b1_max],[b2_min,b2_max],[b3_min,b3_max]]

def generate_grid(b3_max, b1_max=2, b2_max=2):
    l1 = np.linspace(0, b1_max, b1_max+1).tolist()
    l2 = np.linspace(0, b2_max, b2_max+1).tolist()
    l3 = np.linspace(0, b3_max, b3_max+1).tolist()
    return [l1,l2,l3]

def upper_bound(euler):
    return 0.85*euler

def lower_bound(euler):
    return 0.79*euler

def func(p, *d):
    f1, f2, f3 = d
    return (f1(p[0], p[1], p[2]), f2(p[0], p[1], p[2]), f3(p[0], p[1], p[2]))

def constraint(Series, sol, euler, b_list):
    b1, b2, b3 = b_list
    vol = Series(sol[0], sol[1], sol[2])
    vol_min = lower_bound(euler)
    vol_max = upper_bound(euler)
    if 1.0/vol <= vol_max and 1.0/vol >= vol_min:
        return 1, vol
    
    print 'volume: ', vol, ' is out of bounds.'
    
    return 0, -1
    
def NSolve(Series, d, euler, b_list, bound=[[0,1],[0,1],[0,1]], MAX_COUNT=3):
    vol = -1
    sol = -1
    const = 0
    count = 0
    b1, b2, b3 = b_list
    
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
        print 'reset starting point: ', d1_0, d2_0, d3_0
        
        try:
            sol = fsolve(func, x0 = np.array([d1_0, d2_0, d3_0]), args = d)
            print 'solution: ', sol
            print 'guessed vol: ', Series(sol[0], sol[1], sol[2])
        except:
            continue
        
        const, vol = constraint(Series, sol, euler, b_list)

    print 'Done.'

    return vol, sol


def solver(series, h_max, b_list, euler, sol_max=100):
    b1, b2, b3 = b_list
    # find derivative
    d1 = diff(series, b1)
    d2 = diff(series, b2)
    d3 = diff(series, b3)
    d1 = lambdify((b1,b2,b3),d1)
    d2 = lambdify((b1,b2,b3),d2)
    d3 = lambdify((b1,b2,b3),d3)
    d = (d1, d2, d3)
    series = lambdify((b1,b2,b3),series)
    
    # divide solution space into grids
    bounds = generate_grid(h_max)
    for i in range(len(bounds[0])):
        for j in range(len(bounds[1])):
            for k in range(len(bounds[2])):
                bound = find_bound(i,j,k, bounds)
                # try solve
                vol, sol = NSolve(series, d, euler, b_list, bound=bound)
                if type(sol) == int or type(vol) == int:
                    # sol = -1 and vol = -1
                    print ('range ', bound,' does not work')
                    continue
                if type(sol) == np.ndarray:
                    # edge case
                    sol = sol.tolist()
                if sol[0] > sol_max or sol[1] > sol_max or sol[2] > sol_max:
                    print ('solution out of bounds.')
                    continue
                if 0 < vol < 1:
                    print ('vol:', vol, '; sol:', sol)
                    return vol, sol
    print('cannot find solution.')
    return -1,-1

def generate_vol(h_max, coeff=1):
    # Input:
    #   h_max: the max height of the polytope
    #   coeff: an empirical factor that determines the search space
    # Output:
    # triang_list: data for each polytope of height [h1,h2,h3]
    #   [h1,h2,h3]: heights
    #   vol: minimum volume
    #   sol: solution (i.e. the b's in the hilber series)
    #   prism: polytope point set
    #   series: hilbert series
    #   triang: triangulation of the polytope
    #   power: order of the t's in the hilbert series
    triang_list = []
    for h1 in range(1,h_max):
        for h2 in range(h1+1):
            for h3 in range(h2+1):
                print(h1,h2,h3)
                prism, series, triang, power, euler = lift_prism(h1,h2,h3)
                vol, sol = solver(series, h1*coeff, euler)
                triang_list.append([[h1,h2,h3], vol, sol, prism, series, triang, power, euler])
                print('')
    return triang_list

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

def load_series(file_name, home=expanduser("~/Calabi_Yau")):
    path_name = home+'/data/'+file_name
    return load_file(path_name)
    
def find_vol_from_series(file_name, start=0, coeff=1, home=expanduser("~/Calabi_Yau")):
    output_name = 'vol_3913_'+file_name
    output_path = home+'/data/'+output_name
    
    data_list = load_series(file_name)
    b1,b2,b3 = sp.symbols('b1 b2 b3')
    
    counter = 0
    
    for data in data_list:
        counter += 1
        if counter <= start:
            continue
        if counter%10 == 0:
            print counter
        h_list = []
        if type(data[0]) == list:
            h_list = data[0]
        else:
            h_list = eval(data[0])
        print(h_list)
        series = eval(data[2])
        euler = eval(data[-1])
        h1 = max(h_list)
        vol, sol = solver(series, h1*coeff, [b1,b2,b3], euler) #series, h_max, b_list, euler
        write_file(output_path, [h_list, vol, sol])
        print('===========')

        
find_vol_from_series('series_0_49.npy', start=3913) #3913
# print load_series('vol_series_0_2.npy')