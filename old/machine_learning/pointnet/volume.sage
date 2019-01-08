import math
import numpy as np
import scipy
import scipy.optimize
#import matplotlib.pyplot as plt
import csv
#from scipy.optimize import minimize
#import tensorflow as tf

PointConfiguration.set_engine('internal')
#PointConfiguration.set_engine('topcom')
#print 'Done.'

#Helper functions to initialize input data

def reflex_poly(n):

    polytopes_list = [1, 5, 6, 7, 8, 25, 26, 27, 28, 29, 30, 31, 82, 83, 84, 85, 219, 220]

    P = list(ReflexivePolytope(3, polytopes_list[n-1]-1).vertices())
    pts = []

    for j in range(len(P)):
        pts.append(list(P[j]))

    return pts

def reflex_poly_list(num_poly):

    if num_poly > 18:
        print 'Please input integers <= 18'
        return -1

    polytopes_list = [1, 5, 6, 7, 8, 25, 26, 27, 28, 29, 30, 31, 82, 83, 84, 85, 219, 220]

    pts = []

    for i in range(num_poly):
        pts.append(reflex_poly(i))

    print 'Initialize reflexive polytope done.'

    return pts

def input_data(input_path):

    #Code to load 'poly_out.txt'
    with open(input_path) as f:
        pts = []
        pts_tmp = []

        for line in f:
            #print line
            pts_tmp = eval(line)
            pts.append(pts_tmp)

        '''
        #Code to load 'result.txt'
        for line in f:
            pt = line.split(' ')
            print line
            if line == '\n' or pt == '\n':
                pts.append(pts_tmp)
                pts_tmp = []
                continue
            for i in range(len(pt)):
                print pt[i]
                pt[i] = float(pt[i].strip(','))
            pts_tmp.append(pt)
        if line != '\n':
            pts.append(pts_tmp)
        '''

    return pts

#print 'Done.'

#Helper functions to find the Hilbert Series

def exist(pts, latt):
    latt = np.array(latt)
    for i in range(pts.shape[0]):
        if pts[i][0]==latt[0]:
            if pts[i][1]==latt[1]:
                if pts[i][2]==latt[2]:
                    return 1
    return 0

#Compute cross product of three 4-vectors
def four_cross(v1, v2, v3, v4):
    #print "input vectors: ", v1, v2, v3, v4
    v = np.zeros((4,))
    '''
    v1 = list(map(int, v1))
    v2 = list(map(int, v2))
    v3 = list(map(int, v3))
    '''
    counter = 0

    for i in range(4):
        mat = [v1[np.arange(len(v1))!=i].tolist(), v2[np.arange(len(v2))!=i].tolist(), v3[np.arange(len(v3))!=i].tolist()]
        mat = matrix(ZZ, mat)
        #print 'matrix: '
        #print mat
        if counter == 1:
            v[i] = -1*mat.det()
            counter = 0
            #print 'neg: ', v[i]
            continue
        elif counter == 0:
            v[i] = mat.det()
            counter = 1
            #print 'pos: ', v[i]
    #print v
    mat = matrix(RR, [v1.tolist(), v2.tolist(), v3.tolist(), v4.tolist()])

    if mat.det() < 0:
        #print 'original: ', v
        v = -1*v
        #print 'changed: ', v
    #print 'vector: ', v
    return v

'''
def Hilb(tri, p, output):
    num_tri = len(tri)
    len_tri = len(tri[0])
    triang_list = np.zeros((num_tri, len_tri, 4))

    #Convert each element of p into a 4-vector
    #whose last entry equals to 1
    for i in range(num_tri):
        for j in range(len_tri):
            triang_list[i][j] = np.append(np.array(p[tri[i][j]]) , 1)

    #print 'triang_list: '
    #print triang_list   '''
def Hilb(triang_list, output):
    triang = np.array(triang_list)
    power = np.zeros(shape = triang.shape)
    Hilb = 0
    t = var('t')
    t1 = var('t1')
    t2 = var('t2')
    t3 = var('t3')
    t4 = var('t4')
    for tri in range(triang.shape[0]):
        hilb = 1
        t_prod = 1
        for i in range(4):
            #Multiplying by -1 is optional
            power[tri][i] = -1*four_cross(triang[tri][i], triang[tri][np.remainder(i+1, 4)], triang[tri][np.remainder(i+2, 4)], triang[tri][np.remainder(i+3, 4)])
            t_prod = t1**(int(power[tri][i][0]))*t2**(int(power[tri][i][1]))*t3**(int(power[tri][i][2]))*t4**int((power[tri][i][3]))
            hilb *= (1-t_prod)**(-1)
        #print 'Hilbert: ', hilb
        Hilb += hilb
    #print 'Hilb: ', Hilb()
    #print Hilb(t1=t, t2=t, t3=t).series(t4, 3)
    #print "p-q web: ", power


    m = var('m')
    b1 = var('b1')
    b2 = var('b2')
    b3 = var('b3')
    b4 = var('b4')
    Hilb *= m**4

    #print 'Hilb: ', str(Hilb(t1 = (m*b1).exp(), t2 = (m*b2).exp(), t3 = (m*b3).exp(), t4 = (m*4).exp())).replace('e', 'E')


    Series = Hilb(t1 = (m*b1).exp(), t2 = (m*b2).exp(), t3 = (m*b3).exp(), t4 = (m*4).exp()).series(m==0, 1)
    Series = Series.truncate()
    #Series = limit(Hilb(t1 = (m*b1).exp(), t2 = (m*b2).exp(), t3 = (m*b3).exp(), t4 = (m*4).exp()), m=0)
    #print 'Series: ', Series

    if output != 0:
        output.write("%s\n" % Series)

    return Series

def Triang(p, output):
    #Input output = 0 if do not want output
    pts = np.array(p)
    poly = Polyhedron(p)
    pts_max = int(max(np.amax(np.absolute(pts), axis=0)))+1
    pts_new = pts
    for i in range(-pts_max, pts_max):
        for j in range(-pts_max, pts_max):
            for k in range(-pts_max, pts_max):
                latt = [i,j,k]
                if exist(pts, latt)==1:
                    continue
                if poly.contains(latt) == 1:
                    pts_new = np.append(pts_new, np.array(latt).reshape((1,3)), axis = 0)
    #print 'pts_new: ', pts_new
    pts_new = pts_new.tolist()
    points = PointConfiguration(pts_new)
    triang = points.triangulate()
    triang = [list(triang[i]) for i in range(len(triang))]
    triang_new = []
    check_triang(triang, pts_new, triang_new)
    #print 'triangulate: ', triang_new

    #Calculate the Hilbert series
    #Series = Hilb(triang, pts_new, output)
    Series = Hilb(triang_new, output)

    return Series
'''
def check_triang(triang, pts):
    triang_new = triang
    for i in range(len(triang)):
        #triang_new.append(list(triang[i]))
        poly_tmp = Polyhedron(list(triang[i]))
        for j in range(len(pts)):
            if poly_tmp.contains(pts[j]):
                triang_tmp = PointConfiguration(list(triang[i])).triangulate()
                triang_new.remove(triang[i])
                triang_new.append(triang_tmp)
    return triang_new

'''
def check_one_triang(triang, pts):
    triang_list = triang
    #print 'pts: ', pts
    #print 'triang: ', triang
    vert = [pts[k] for k in triang]
    poly = Polyhedron(vert)

    for h in range(len(pts)):
        if h in triang:
            continue
        elif poly.contains(pts[h]):
            triang_list.append(h)

    return triang_list

def check_triang(triang, pts, triang_new):
    #triang_new = []
    for i in range(len(triang)):
        check = check_one_triang(triang[i], pts)
        #print "check ", i, ": ", check
        len_check = len(check)

        if len_check <= 4:
            triang_save = [np.append(np.array(pts[k]),1).tolist() for k in check]
            #print 'add: ', triang_save
            triang_new.append(triang_save)

        else:
            pts_tmp = [pts[j] for j in check]
            #print "vertices: ", pts_tmp
            triang_tmp = PointConfiguration(pts_tmp).triangulate()
            triang_tmp = [list(triang_tmp[i]) for i in range(len(triang_tmp))]
            #print "further triang: ", triang_tmp
            check_triang(triang_tmp, pts_tmp, triang_new)
    #print triang_new


#print 'Done.'

#Helper functions to find the minimum of the volume function
#(None of these work as well as Mathematica yet.)

from scipy.optimize import fsolve
from sympy import nsolve

def min_function(Series, p):
    try:
        return (Series(b1=p[0], b2=p[1], b3=p[2]))^2
    except:
        return -1

def min_constraint(Series, a, b, c):
    if min_function(Series, [a,b,c])== -1:
        print 'bad try: ', a,' ', b , ' ', c
        return -100
    try:
        print 'val: ', Series(b1=a, b2=b, b3=c)
        return 1-Series(b1=a, b2=b, b3=c)
    except:
        return -100

def Find_Minimum(Series):
    function = lambda p: min_function(Series, p)

    constraint = ({'type': 'ineq', 'fun': lambda p:  min_constraint(Series, p[0], p[1], p[2])})
    #constraint = ({'type': 'ineq', 'fun': lambda p:  1-Series(b1=p[0], b2=p[1], b3=p[2])})

    solution = scipy.optimize.minimize(function, (0.7,0.2,0.3), constraints=constraint)

    return solution

def eval_sol(Series, sol):
    try:
        return Series(b1 = sol[0].rhs(), b2 = sol[1].rhs(), b3 = sol[2].rhs())
    except:
        print 'no sol'
        return CDF(I)

def D_Minimum(Series):
    d1 = diff(Series, b1)
    d2 = diff(Series, b2)
    d3 = diff(Series, b3)

    solution = solve([d1 == 0, d2 == 0, d3 == 0], b1, b2, b3)
    #print 'solution: ', solution

    sol_len = len(solution)

    new_sol = []

    for j in range(sol_len):
        if solution[j][0].rhs() in RR and solution[j][1].rhs() in RR and solution[j][2].rhs() in RR:
            vol = eval_sol(Series, solution[j])
            print vol
            if vol not in RR:
                continue
            new_sol.append(vol)
            #new_sol.append(solution[j])
    #solution = new_sol

    print 'new_sol: ', new_sol

    vol_min_abs = min(new_sol)

    return vol_min_abs

def func(p, *d):
    f1, f2, f3 = d
    return (f1(b1 = p[0], b2 = p[1], b3 = p[2]), f2(b1 = p[0], b2 = p[1], b3 = p[2]), f3(b1 = p[0], b2 = p[1], b3 = p[2]))

def constraint(Series, sol):
    vol = Series(b1 = sol[0], b2 = sol[1], b3 = sol[2])
    if vol <= 1 and vol >= 0:
        return 1, vol

    print 'volume: ', vol, ' is out of bounds.'

    return 0, -1

def NSolve(Series):
    d1 = diff(Series, b1)
    d2 = diff(Series, b2)
    d3 = diff(Series, b3)
    d = (d1, d2, d3)
    const = 0
    count = 0

    while const == 0:
        d1_0 = np.random.uniform(low=-10, high=10)
        d2_0 = np.random.uniform(low=-10, high=10)
        d3_0 = np.random.uniform(low=-10, high=10)
        print 'reset starting point: ', d1_0, d2_0, d3_0
        try:
            sol = fsolve(func, x0 = np.array([d1_0, d2_0, d3_0]), args = d)
            print 'solution: ', sol
        except:
            continue
        if abs(sol[0]) > 100 or abs(sol[1]) > 100 or abs(sol[2]) > 100:
            continue
        const, vol = constraint(Series, sol)

        count += 1
        if count > 1000:
            print 'Infinite loop. Force stop.'
            return -1, -1

    print 'Done.'

    return vol, sol

#print 'Done.'

#Generate hilbert series and volume

MAX_PTS = 35

def generate_hilbert_vol(input_path, output_hilb_path, output_vol_path):
    pts = input_data(input_path)
    for i in range(len(pts)):
        pts_new = pts[i]
        points = PointConfiguration(pts_new)

        #Triangulate
        triang = points.triangulate()
        triang = [list(triang[i]) for i in range(len(triang))]
        triang_new = []
        check_triang(triang, pts_new, triang_new)

        #Calculate Hilbert Series (write to output)
        Series = Hilb(triang_new, 0)

        #Calculate Volume (write to output)
        vol, sol = NSolve(Series)

        num_pts = len(pts_new)

        if num_pts < MAX_PTS:
            for j in range(MAX_PTS - num_pts):
                pts_new.append([1j,1j,1j])

        output_list = [pts_new, vol]
        print i,'-th output: ', output_list
        output_hilb = open(output_hilb_path, 'a')
        output_vol = open(output_vol_path, 'a')
        output_hilb.write("%s\n" % str([Series, sol]))
        output_vol.write("%s\n" % str(output_list))
        output_hilb.close()
        output_vol.close()

#output/polygon/poly_out_2.txt
#output/series/series_cube_30.txt
#output/train/train.txt

input_path = raw_input("Please input the input path (example: output/polygon/poly_out_1024.txt): ")
input_path = str(input_path)
output_hilb_path = raw_input("Please input the output path of the hilbert series (example: output/series/series_cube_1024.txt): ")
output_hilb_path = str(output_hilb_path)
output_vol_path = raw_input("Please input the output path of the volume (example: output/vol/vol_cube_30.txt): ")
output_vol_path = str(output_vol_path)

#output_hilb_path = 'output/series/output_2.txt'
#output_vol_path = 'output/vol/output_cube_30.txt'
#output_hilb = open(output_hilb_path, 'w')
#output_vol = open(output_vol_path, 'w')
generate_hilbert_vol(input_path, output_hilb_path, output_vol_path)
#output_hilb.close()
#output_vol.close()
print 'Done.'
