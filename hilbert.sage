import math
import numpy as np
import scipy
#import matplotlib.pyplot as plt
import csv

PointConfiguration.set_engine('internal')

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
    v = np.zeros((4,))
    counter = 0
    #print 'input vector: ', v1, v2, v3
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
    #print triang_list
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
            t_prod = t1^(power[tri][i][0])*t2^(power[tri][i][1])*t3^(power[tri][i][2])*t4^(power[tri][i][3])
            hilb *= (1-t_prod)^(-1)
        #print 'Hilbert: ', hilb
        Hilb += hilb
    #print 'Hilb: ', Hilb(t1 = t, t2 = t, t3 = t, t4 = t^4).factor()
    #print Hilb(t1=t, t2=t, t3=t).series(t4, 3)
    output.write("p-q web: %s\n" % power)
    
    m = var('m')
    b1 = var('b1')
    b2 = var('b2')
    b3 = var('b3')
    b4 = var('b4')
    Hilb *= m^4
    
    Series = Hilb(t1 = (m*b1).exp(), t2 = (m*b2).exp(), t3 = (m*b3).exp(), t4 = (m*4).exp()).series(m==0, 1)
    Series = Series.truncate()
    d1 = diff(Series, b1)
    d2 = diff(Series, b2)
    d3 = diff(Series, b3)
    
    solution = solve([d1 == 0, d2 == 0, d3 == 0], b1, b2, b3)
    
    sol_len = len(solution)
    
    vol_min_abs = Series(b1 = solution[0][0].rhs(), b2 = solution[0][1].rhs(), b3 = solution[0][2].rhs())
    
    if sol_len > 1:
        #print 'sol_len: ', sol_len
        for i in range(1, sol_len):
            vol_min = Series(b1 = solution[i][0].rhs(), b2 = solution[i][1].rhs(), b3 = solution[i][2].rhs())
            if vol_min < vol_min_abs and vol_min > 0:
                vol_min_abs = vol_min
    
    return vol_min_abs
    
def Triang(p, output):
    pts = np.array(p)
    poly = Polyhedron(p)
    pts_max = int(max(np.amax(pts, axis=0)))
    pts_new = pts
    for i in range(0, pts_max):
        for j in range(0, pts_max):
            for k in range(0, pts_max):
                latt = [i,j,k]
                if exist(pts, latt)==1:
                    continue
                if poly.contains(latt) == 1:
                    pts_new = np.append(pts_new, np.array(latt).reshape((1,3)), axis = 0)  
    #print 'pts_new: ', pts_new
    pts_new = pts_new.tolist()
    points = PointConfiguration(pts_new)
    triang = points.triangulate()
    triang = list(triang)
    #print 'triangulate: ', triang
    
    #Calculate the Hilbert series
    vol_min = Hilb(triang, pts_new, output)
    output.write("Vol: %s\n" % vol_min)
    #print 'vol_min: ', vol_min

input_path = raw_input("Please input the input path: ")
input_path = str(input_path)
output_path = raw_input("Please input the output path: ")
output_path = str(output_path)

#input_path = 'input.csv'
#output_path = 'output.txt'
with open(input_path) as f:
    pts = []
    pts_tmp = []
    
    for line in f:
        pt = line.split(' ')
        if line == '\n' or pt == '\n':
            pts.append(pts_tmp)
            pts_tmp = []
            continue
        for i in range(len(pt)):
            pt[i] = float(pt[i].strip(','))
        pts_tmp.append(pt)
    if line != '\n':
        pts.append(pts_tmp)
        
output = open(output_path, 'w')

for idx in range(len(pts)):
    output.write("polytope %s\n" % idx)
    Triang(pts[idx], output)
    output.write("\n")

output.close()
