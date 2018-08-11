import math
import numpy as np
import scipy
import scipy.optimize
import csv
from scipy.optimize import fsolve

PointConfiguration.set_engine('internal')
#PointConfiguration.set_engine('topcom')
print 'Done.'

def exist(pts, latt):
    latt = np.array(latt)
    for i in range(pts.shape[0]):
        if pts[i][0]==latt[0]:
            if pts[i][1]==latt[1]:
                if pts[i][2]==latt[2]:
                    return 1
    return 0

def contain(poly, latt):
    if poly.contains(latt) == 1:
        return 1
    else:
        poly_latt = Polyhedron(vertices = [tuple(latt)])
        vert = next(poly_latt.vertex_generator())
        face_eq = poly.Hrepresentation()
        for eq in face_eq:
            if eq.contains(vert) != 1:
                return 0
        return 1
    return 0

def check_latt(p):
    pts = np.array(p)
    pts_max = int(max(np.amax(pts, axis=0)))+1
    pts_min = int(min(np.amin(pts, axis=0)))-1
    print 'pts_max and pts_min: ', pts_max, pts_min
    poly = Polyhedron(p)
    pts_new = pts
    for i in range(pts_min, pts_max):
        for j in range(pts_min, pts_max):
            for k in range(pts_min, pts_max):
                latt = [i,j,k]
                if exist(pts, latt) == 1:
                    continue
                if contain(poly, latt) == 1:
                    pts_new = np.append(pts_new, np.array(latt).reshape((1,3)), axis = 0)  
    pts_new = pts_new.tolist()
    return pts_new


def four_cross(v1, v2, v3, v4):
    #Compute cross product of three 4-vectors
    #print "input vectors: ", v1, v2, v3, v4
    v = np.zeros((4,))
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

def Hilb(triang_list):
    triang = []
    # Add 1 at the end of all verticies
    for tetra in triang_list:
        tetra_new = []
        for vert in tetra:
            vert_new = np.append(np.array(vert),1).tolist()
            tetra_new.append(vert_new)
        triang.append(tetra_new)
    triang = np.array(triang)
    #print 'input: ', triang_list
    #print 'triang: ', triang
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
            t_prod = t1^(int(power[tri][i][0]))*t2^(int(power[tri][i][1]))*t3^(int(power[tri][i][2]))*t4^int((power[tri][i][3]))
            hilb *= (1-t_prod)^(-1)
        #print 'Hilbert: ', hilb
        Hilb += hilb
    
    
    m = var('m')
    b1 = var('b1')
    b2 = var('b2')
    b3 = var('b3')
    b4 = var('b4')
    Hilb *= m^4
    
    #print 'Hilb: ', str(Hilb(t1 = (m*b1).exp(), t2 = (m*b2).exp(), t3 = (m*b3).exp(), t4 = (m*(2-(b1+b2+b3))).exp())).replace('e', 'E')
    
    
    Series = Hilb(t1 = (m*b1).exp(), t2 = (m*b2).exp(), t3 = (m*b3).exp(), t4 = (m*(2-(b1+b2+b3))).exp()).series(m==0, 1)
    Series = Series.truncate()
    
    return Series
    
# Triangulate Triangular Prism
def init_prism(h):
    if h == 0:
        prism = []
        series = 0
        return prism, series
    prism_1 = [[[0,0,0],[1,0,0],[0,1,0],[1,0,1]]]
    prism_2 = [[[0,0,0],[1,0,1],[0,1,1],[0,1,0]]]
    prism_3 = [[[0,0,0],[1,0,1],[0,0,1],[0,1,1]]]
    prism = prism_1 + prism_2 + prism_3
    if h == 1:
        series = Hilb(prism)
        return prism, series
    prism_stack = prism
    for z in range(1, h):
        move = np.array([0,0,z])
        for i in range(3):
            # For i-th tetrahedron
            tetra_moved = []
            for j in range(4):
                # For j-th vertex:
                vert_moved = (np.array(prism[i][j])+move).tolist()
                tetra_moved.append(vert_moved)
            prism_stack.append(tetra_moved)
    # Also initialize the Hilbert series
    series = Hilb(prism_stack)
    return prism_stack, series

def Triang_prism(h1, h2, h3):
    # Phase 1: Cube blocks up to min(h1, h2, h3)
    h_list = [h1, h2, h3]
    h_min = min(h_list)
    h_list.remove(h_min)
    prism, series = init_prism(h_min)
    # Phase 2: Tetra with min pt as its apex
    h_min_2 = min(h_list)
    for i in range(h_min_2 - h_min):
        # Also find the Hilbert series for each tetrahedron
        prism_1 = [[[1,0,h_min],[0,0,h_min+i],[0,1,h_min+i],[0,1,h_min+i+1]]]
        series_1 = Hilb(prism_1)
        prism_2 = [[[1,0,h_min],[0,0,h_min+i],[0,0,h_min+i+1],[0,1,h_min+i+1]]]
        series_2 = Hilb(prism_2)
        prism += prism_1 + prism_2
        series += series_1 + series_2
    
    # Phase 3: Tetra with max pt as its apex
    h_max = max(h_list)
    for i in range(h_max-h_min_2):
        prism_new = [[[1,0,h_min],[0,0,h_min_2+i],[0,0,h_min_2+i+1],[0,1,h_min_2]]]
        series_new = Hilb(prism_new)
        prism += prism_new
        series += series_new

    return prism, series
    
def scipy_opt(series, num_sec, SIDE_LENGTH, out_file):
    d1 = diff(series, b1)
    d2 = diff(series, b2)
    d3 = diff(series, b3)
    print 'd1: ', d1
    fun = lambda b: (d1(b1=b[0], b2=b[1], b3=b[2]))**2 + (d2(b1=b[0], b2=b[1], b3=b[2]))**2 + (d3(b1=b[0], b2=b[1], b3=b[2]))**2
    bnds = ((0,1),(0,1),(0,1))
    cons = [{'type': 'ineq', 'fun': lambda b:  1-(2-(b[0]+b[1]+b[2]))},
            {'type': 'ineq', 'fun': lambda b:  2-(b[0]+b[1]+b[2])}]
    
    d1_min = 0
    d1_max = 1
    d2_min = 0
    d2_max = 1
    d3_min = 0
    d3_max = 1
    d1_0 = np.random.uniform(low=d1_min, high=d1_max)
    d2_0 = np.random.uniform(low=d2_min, high=d2_max)
    d3_0 = np.random.uniform(low=d3_min, high=d3_max)
    b_init = [d1_0, d2_0, d3_0]
    print 'initial guess: ', b_init
    res = scipy.optimize.minimize(fun, b_init, method='SLSQP', bounds=bnds, constraints=cons)
    return res

def func(p, *d):
    f1, f2, f3 = d
    return (f1(b1 = p[0], b2 = p[1], b3 = p[2]), f2(b1 = p[0], b2 = p[1], b3 = p[2]), f3(b1 = p[0], b2 = p[1], b3 = p[2]))

def constraint(Series, sol, SIDE_LENGTH, out_file):
    vol = -1*Series(b1 = sol[0], b2 = sol[1], b3 = sol[2])
    if vol > 0 and vol <= (3*SIDE_LENGTH)**3:
        return 1, vol
    
    out_file.write('volume: %f is out of bounds.' % (vol))
    print 'volume: ', vol, ' is out of bounds.'
    
    return 0, -1

def constraint_tmp(Series, sol, SIDE_LENGTH, out_file):
    vol = Series(b1 = sol[0], b2 = sol[1], b3 = sol[2])
    if vol <= 1 and vol >= 1/((3*SIDE_LENGTH)**3):
        return 1, vol

    print 'volume: ', vol, ' is out of bounds.'

    return 0, -1

def NSolve(Series, d, SIDE_LENGTH, bound, out_file):
    const = 0
    count = 0
    MAX_COUNT = 2
    
    d1_min = bound[0][0]
    d1_max = bound[0][1]
    d2_min = bound[1][0]
    d2_max = bound[1][1]
    d3_min = bound[2][0]
    d3_max = bound[2][1]
    
    while const == 0:
        if count >= MAX_COUNT:
            print 'No good'
            return -1,-1
        d1_0 = np.random.uniform(low=d1_min, high=d1_max)
        d2_0 = np.random.uniform(low=d2_min, high=d2_max)
        d3_0 = np.random.uniform(low=d3_min, high=d3_max)
        print 'reset starting point: ', d1_0, d2_0, d3_0
        count += 1

        #try:
        sol = fsolve(func, x0 = np.array([d1_0, d2_0, d3_0]), args = d)
        if type(sol) == np.ndarray:
            sol = sol.tolist()
        print 'solution: ', sol
        b4 = 2 - (sol[0] + sol[1] + sol[2]) 
        if sol[0] > 1 or sol[1] > 1 or sol[2] > 1 or b4 > 1 or sol[0] < 0 or sol[1] < 0 or sol[2] < 0 or b4 < 0:
            print 'solution ', [sol, b4], ' is out of bounds.'
            print 'wrong vol: ', Series(b1 = sol[0], b2 = sol[1], b3 = sol[2])
            continue
        #except:
        #    continue
        
        const, vol = constraint(Series, sol, SIDE_LENGTH, out_file)
        

    print 'Done.'

    return vol, sol

def grid_NSolve(Series, d, num_sec, SIDE_LENGTH, out_file):
    length = 1/float(num_sec)
    counter = 0
    vol_list = []
    sol_list = []
    const = 0
    count = 0
    
    for i in range(num_sec):
        for j in range(num_sec):
            for k in range(num_sec):
                print 'grid number: ', counter+1
                out_file.write('grid #%d\n' % (counter+1))
                counter += 1
                b1_range = [i*length,(i+1)*length]
                b2_range = [j*length,(j+1)*length]
                b3_range = [k*length,(k+1)*length]
                bound = [b1_range, b2_range, b3_range]
                print 'bounds: '
                print bound
                out_file.write('grid bound: %s\n' % (str(bound)))
                vol,sol = NSolve(Series, d, SIDE_LENGTH, bound, out_file)
                try:
                    if sol == -1:
                        print 'grid ', counter, ' dose not work.'
                        continue
                except:
                    if vol != -1:
                        print 'result vol: ', vol
                        print 'result sol: ', sol
                        out_file.write('vol: %f\nsol: %s\n' % (vol, str(sol)))
                        if type(sol) == np.ndarray:
                            sol = sol.tolist()
                        vol_list.append(vol)
                        sol_list.append(sol)
                        #return vol, sol
                
    print 'num blocks: ', counter
    
    return vol_list, sol_list

def iter_grid_NSolve(Series, max_sec, SIDE_LENGTH, out_file):
    d1 = diff(Series, b1)
    d2 = diff(Series, b2)
    d3 = diff(Series, b3)
    d = (d1, d2, d3)
    vol_min_list = []
    sol_min_list = []
    for i in range(3, max_sec+1):
        print 'Try grid with ', i, ' length segments.'
        out_file.write('grid with length #%d\n' % i)
        vol_list, sol_list = grid_NSolve(Series, d, i, SIDE_LENGTH, out_file)
        if len(vol_list) != 0 and len(sol_list) != 0:
            print 'vol_list ', i, ': ', vol_list
            vol_min = min(vol_list)
            sol_min = sol_list[np.argmin(vol_list)]
            if type(sol_min) == np.ndarray:
                sol_min = sol_min.tolist()
            print 'vol min: ', vol_min
            vol_min_list.append(vol_min)
            sol_min_list.append(sol_min)
    if len(vol_min_list) != 0 and len(sol_min_list) != 0:
        print 'vol_min_list: ', vol_min_list
        print 'sol_min_list: ', sol_min_list
        vol_min = min(vol_min_list)
        sol_min = sol_min_list[np.argmin(vol_min_list)]
        return vol_min, sol_min
    else:
        return -1, -1

SIDE_LENGTH = 5
b1 = var('b1')
b2 = var('b2')
b3 = var('b3')
b4 = var('b4')
'''
series = 1/((8*b1 + 5*b2 + 6*b3 - 10)*(7*b1 + 6*b2 + 6*b3 - 10)*(b1 + 2*b2 + b3 - 2)*b1) + 1/((6*b1 + 4*b2 + 5*b3 - 8)*(5*b1 + 5*b2 + 5*b3 - 8)*(b1 + 2*b2 + b3 - 2)*b1) + 1/((4*b1 + 3*b2 + 4*b3 - 6)*(3*b1 + 4*b2 + 4*b3 - 6)*(b1 + 2*b2 + b3 - 2)*b1) + 1/((2*b1 + 2*b2 + 3*b3 - 4)*(b1 + 3*b2 + 3*b3 - 4)*(b1 + 2*b2 + b3 - 2)*b1) + 1/((7*b1 + 6*b2 + 6*b3 - 10)*(6*b1 + 4*b2 + 5*b3 - 8)*(b1 - b2)*b1) + 1/((5*b1 + 5*b2 + 5*b3 - 8)*(4*b1 + 3*b2 + 4*b3 - 6)*(b1 - b2)*b1) + 1/((3*b1 + 4*b2 + 4*b3 - 6)*(2*b1 + 2*b2 + 3*b3 - 4)*(b1 - b2)*b1) + 1/((2*b1 - b3)*(b1 - b2)*(b1 - 2*b2 - 2*b3 + 2)*b1) + 1/((b1 + 3*b2 + 3*b3 - 4)*(b1 - b2)*b1*(b2 + 2*b3 - 2)) - 1/((b1 + 2*b2 + b3 - 2)*(b1 - 2*b2 - 2*b3 + 2)*b1*(b2 + 2*b3 - 2)) - 1/((10*b1 + 3*b2 + 6*b3 - 10)*(8*b1 + 2*b2 + 5*b3 - 8)*(b1 - b2)*b2) - 1/((8*b1 + 2*b2 + 5*b3 - 8)*(6*b1 + b2 + 4*b3 - 6)*(b1 - b2)*b2) - 1/((6*b1 + b2 + 4*b3 - 6)*(4*b1 + 3*b3 - 4)*(b1 - b2)*b2) - 1/((4*b1 + 3*b3 - 4)*(2*b1 - b2 + 2*b3 - 2)*(b1 - b2)*b2) + 1/((2*b1 - b2 + 2*b3 - 2)*(b1 - b2)*(2*b2 - b3)*b2) - 1/((4*b1 + 4*b2 + b3 - 4)*(2*b1 - b3)*(2*b2 - b3)*b3)
'''
series = -1/2/((4*b1 + 4*b2 + 3*b3 - 4)*(b1 + b2 + b3 - 1)*b1*b2) - 1/2/((b1 + b2 + b3 - 1)*b1*b2*b3)

h1 = 5
h2 = 5
h3 = 2
prism, series = Triang_prism(h1,h2,h3)
print 'series for: ', h1, h2, h3
print series
out_path = '%d_%d_%d_test.txt' % (h1,h2,h3)
out_file = open(out_path, 'w')
vol, sol = iter_grid_NSolve(series, 4, SIDE_LENGTH, out_file)
out_file.close()