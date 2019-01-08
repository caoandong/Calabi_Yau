
# This file was *autogenerated* from the file numerical.sage
from sage.all_cmdline import *   # import sage library

_sage_const_3 = Integer(3); _sage_const_2 = Integer(2); _sage_const_1 = Integer(1); _sage_const_0 = Integer(0); _sage_const_5 = Integer(5); _sage_const_4 = Integer(4)
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
    for i in range(pts.shape[_sage_const_0 ]):
        if pts[i][_sage_const_0 ]==latt[_sage_const_0 ]:
            if pts[i][_sage_const_1 ]==latt[_sage_const_1 ]:
                if pts[i][_sage_const_2 ]==latt[_sage_const_2 ]:
                    return _sage_const_1 
    return _sage_const_0 

def contain(poly, latt):
    if poly.contains(latt) == _sage_const_1 :
        return _sage_const_1 
    else:
        poly_latt = Polyhedron(vertices = [tuple(latt)])
        vert = next(poly_latt.vertex_generator())
        face_eq = poly.Hrepresentation()
        for eq in face_eq:
            if eq.contains(vert) != _sage_const_1 :
                return _sage_const_0 
        return _sage_const_1 
    return _sage_const_0 

def check_latt(p):
    pts = np.array(p)
    pts_max = int(max(np.amax(pts, axis=_sage_const_0 )))+_sage_const_1 
    pts_min = int(min(np.amin(pts, axis=_sage_const_0 )))-_sage_const_1 
    print 'pts_max and pts_min: ', pts_max, pts_min
    poly = Polyhedron(p)
    pts_new = pts
    for i in range(pts_min, pts_max):
        for j in range(pts_min, pts_max):
            for k in range(pts_min, pts_max):
                latt = [i,j,k]
                if exist(pts, latt) == _sage_const_1 :
                    continue
                if contain(poly, latt) == _sage_const_1 :
                    pts_new = np.append(pts_new, np.array(latt).reshape((_sage_const_1 ,_sage_const_3 )), axis = _sage_const_0 )  
    pts_new = pts_new.tolist()
    return pts_new


def four_cross(v1, v2, v3, v4):
    #Compute cross product of three 4-vectors
    #print "input vectors: ", v1, v2, v3, v4
    v = np.zeros((_sage_const_4 ,))
    counter = _sage_const_0 
    
    for i in range(_sage_const_4 ):
        mat = [v1[np.arange(len(v1))!=i].tolist(), v2[np.arange(len(v2))!=i].tolist(), v3[np.arange(len(v3))!=i].tolist()]
        mat = matrix(ZZ, mat)
        #print 'matrix: '
        #print mat
        if counter == _sage_const_1 :
            v[i] = -_sage_const_1 *mat.det()
            counter = _sage_const_0 
            #print 'neg: ', v[i]
            continue
        elif counter == _sage_const_0 :
            v[i] = mat.det()
            counter = _sage_const_1 
            #print 'pos: ', v[i]
    #print v
    mat = matrix(RR, [v1.tolist(), v2.tolist(), v3.tolist(), v4.tolist()])
    
    if mat.det() < _sage_const_0 :
        #print 'original: ', v
        v = -_sage_const_1 *v
        #print 'changed: ', v
    #print 'vector: ', v
    return v

def Hilb(triang_list):
    triang = []
    # Add 1 at the end of all verticies
    for tetra in triang_list:
        tetra_new = []
        for vert in tetra:
            vert_new = np.append(np.array(vert),_sage_const_1 ).tolist()
            tetra_new.append(vert_new)
        triang.append(tetra_new)
    triang = np.array(triang)
    #print 'input: ', triang_list
    #print 'triang: ', triang
    power = np.zeros(shape = triang.shape)
    Hilb = _sage_const_0 
    t = var('t')
    t1 = var('t1')
    t2 = var('t2')
    t3 = var('t3')
    t4 = var('t4')
    for tri in range(triang.shape[_sage_const_0 ]):
        hilb = _sage_const_1 
        t_prod = _sage_const_1 
        for i in range(_sage_const_4 ):
            #Multiplying by -1 is optional
            power[tri][i] = -_sage_const_1 *four_cross(triang[tri][i], triang[tri][np.remainder(i+_sage_const_1 , _sage_const_4 )], triang[tri][np.remainder(i+_sage_const_2 , _sage_const_4 )], triang[tri][np.remainder(i+_sage_const_3 , _sage_const_4 )])
            t_prod = t1**(int(power[tri][i][_sage_const_0 ]))*t2**(int(power[tri][i][_sage_const_1 ]))*t3**(int(power[tri][i][_sage_const_2 ]))*t4**int((power[tri][i][_sage_const_3 ]))
            hilb *= (_sage_const_1 -t_prod)**(-_sage_const_1 )
        #print 'Hilbert: ', hilb
        Hilb += hilb
    
    
    m = var('m')
    b1 = var('b1')
    b2 = var('b2')
    b3 = var('b3')
    b4 = var('b4')
    Hilb *= m**_sage_const_4 
    
    #print 'Hilb: ', str(Hilb(t1 = (m*b1).exp(), t2 = (m*b2).exp(), t3 = (m*b3).exp(), t4 = (m*(2-(b1+b2+b3))).exp())).replace('e', 'E')
    
    
    Series = Hilb(t1 = (m*b1).exp(), t2 = (m*b2).exp(), t3 = (m*b3).exp(), t4 = (m*(_sage_const_2 -(b1+b2+b3))).exp()).series(m==_sage_const_0 , _sage_const_1 )
    Series = Series.truncate()
    
    return Series
    
# Triangulate Triangular Prism
def init_prism(h):
    if h == _sage_const_0 :
        prism = []
        series = _sage_const_0 
        return prism, series
    prism_1 = [[[_sage_const_0 ,_sage_const_0 ,_sage_const_0 ],[_sage_const_1 ,_sage_const_0 ,_sage_const_0 ],[_sage_const_0 ,_sage_const_1 ,_sage_const_0 ],[_sage_const_1 ,_sage_const_0 ,_sage_const_1 ]]]
    prism_2 = [[[_sage_const_0 ,_sage_const_0 ,_sage_const_0 ],[_sage_const_1 ,_sage_const_0 ,_sage_const_1 ],[_sage_const_0 ,_sage_const_1 ,_sage_const_1 ],[_sage_const_0 ,_sage_const_1 ,_sage_const_0 ]]]
    prism_3 = [[[_sage_const_0 ,_sage_const_0 ,_sage_const_0 ],[_sage_const_1 ,_sage_const_0 ,_sage_const_1 ],[_sage_const_0 ,_sage_const_0 ,_sage_const_1 ],[_sage_const_0 ,_sage_const_1 ,_sage_const_1 ]]]
    prism = prism_1 + prism_2 + prism_3
    if h == _sage_const_1 :
        series = Hilb(prism)
        return prism, series
    prism_stack = prism
    for z in range(_sage_const_1 , h):
        move = np.array([_sage_const_0 ,_sage_const_0 ,z])
        for i in range(_sage_const_3 ):
            # For i-th tetrahedron
            tetra_moved = []
            for j in range(_sage_const_4 ):
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
        prism_1 = [[[_sage_const_1 ,_sage_const_0 ,h_min],[_sage_const_0 ,_sage_const_0 ,h_min+i],[_sage_const_0 ,_sage_const_1 ,h_min+i],[_sage_const_0 ,_sage_const_1 ,h_min+i+_sage_const_1 ]]]
        series_1 = Hilb(prism_1)
        prism_2 = [[[_sage_const_1 ,_sage_const_0 ,h_min],[_sage_const_0 ,_sage_const_0 ,h_min+i],[_sage_const_0 ,_sage_const_0 ,h_min+i+_sage_const_1 ],[_sage_const_0 ,_sage_const_1 ,h_min+i+_sage_const_1 ]]]
        series_2 = Hilb(prism_2)
        prism += prism_1 + prism_2
        series += series_1 + series_2
    
    # Phase 3: Tetra with max pt as its apex
    h_max = max(h_list)
    for i in range(h_max-h_min_2):
        prism_new = [[[_sage_const_1 ,_sage_const_0 ,h_min],[_sage_const_0 ,_sage_const_0 ,h_min_2+i],[_sage_const_0 ,_sage_const_0 ,h_min_2+i+_sage_const_1 ],[_sage_const_0 ,_sage_const_1 ,h_min_2]]]
        series_new = Hilb(prism_new)
        prism += prism_new
        series += series_new

    return prism, series
    
def scipy_opt(series, num_sec, SIDE_LENGTH, out_file):
    d1 = diff(series, b1)
    d2 = diff(series, b2)
    d3 = diff(series, b3)
    print 'd1: ', d1
    fun = lambda b: (d1(b1=b[_sage_const_0 ], b2=b[_sage_const_1 ], b3=b[_sage_const_2 ]))**_sage_const_2  + (d2(b1=b[_sage_const_0 ], b2=b[_sage_const_1 ], b3=b[_sage_const_2 ]))**_sage_const_2  + (d3(b1=b[_sage_const_0 ], b2=b[_sage_const_1 ], b3=b[_sage_const_2 ]))**_sage_const_2 
    bnds = ((_sage_const_0 ,_sage_const_1 ),(_sage_const_0 ,_sage_const_1 ),(_sage_const_0 ,_sage_const_1 ))
    cons = [{'type': 'ineq', 'fun': lambda b:  _sage_const_1 -(_sage_const_2 -(b[_sage_const_0 ]+b[_sage_const_1 ]+b[_sage_const_2 ]))},
            {'type': 'ineq', 'fun': lambda b:  _sage_const_2 -(b[_sage_const_0 ]+b[_sage_const_1 ]+b[_sage_const_2 ])}]
    
    d1_min = _sage_const_0 
    d1_max = _sage_const_1 
    d2_min = _sage_const_0 
    d2_max = _sage_const_1 
    d3_min = _sage_const_0 
    d3_max = _sage_const_1 
    d1_0 = np.random.uniform(low=d1_min, high=d1_max)
    d2_0 = np.random.uniform(low=d2_min, high=d2_max)
    d3_0 = np.random.uniform(low=d3_min, high=d3_max)
    b_init = [d1_0, d2_0, d3_0]
    print 'initial guess: ', b_init
    res = scipy.optimize.minimize(fun, b_init, method='SLSQP', bounds=bnds, constraints=cons)
    return res

def func(p, *d):
    f1, f2, f3 = d
    return (f1(b1 = p[_sage_const_0 ], b2 = p[_sage_const_1 ], b3 = p[_sage_const_2 ]), f2(b1 = p[_sage_const_0 ], b2 = p[_sage_const_1 ], b3 = p[_sage_const_2 ]), f3(b1 = p[_sage_const_0 ], b2 = p[_sage_const_1 ], b3 = p[_sage_const_2 ]))

def constraint(Series, sol, SIDE_LENGTH, out_file):
    vol = -_sage_const_1 *Series(b1 = sol[_sage_const_0 ], b2 = sol[_sage_const_1 ], b3 = sol[_sage_const_2 ])
    if vol > -_sage_const_1  and vol <= (_sage_const_3 *SIDE_LENGTH)**_sage_const_3 :
        return _sage_const_1 , vol
    
    out_file.write('volume: %f is out of bounds.' % (vol))
    print 'volume: ', vol, ' is out of bounds.'
    
    return _sage_const_0 , -_sage_const_1 

def constraint_tmp(Series, sol, SIDE_LENGTH, out_file):
    vol = Series(b1 = sol[_sage_const_0 ], b2 = sol[_sage_const_1 ], b3 = sol[_sage_const_2 ])
    if vol <= _sage_const_1  and vol >= _sage_const_1 /((_sage_const_3 *SIDE_LENGTH)**_sage_const_3 ):
        return _sage_const_1 , vol

    print 'volume: ', vol, ' is out of bounds.'

    return _sage_const_0 , -_sage_const_1 

def NSolve(Series, d, SIDE_LENGTH, bound, out_file):
    const = _sage_const_0 
    count = _sage_const_0 
    MAX_COUNT = _sage_const_2 
    
    d1_min = bound[_sage_const_0 ][_sage_const_0 ]
    d1_max = bound[_sage_const_0 ][_sage_const_1 ]
    d2_min = bound[_sage_const_1 ][_sage_const_0 ]
    d2_max = bound[_sage_const_1 ][_sage_const_1 ]
    d3_min = bound[_sage_const_2 ][_sage_const_0 ]
    d3_max = bound[_sage_const_2 ][_sage_const_1 ]
    
    while const == _sage_const_0 :
        if count >= MAX_COUNT:
            print 'No good'
            return -_sage_const_1 ,-_sage_const_1 
    
        d1_0 = np.random.uniform(low=d1_min, high=d1_max)
        d2_0 = np.random.uniform(low=d2_min, high=d2_max)
        d3_0 = np.random.uniform(low=d3_min, high=d3_max)
        print 'reset starting point: ', d1_0, d2_0, d3_0
        count += _sage_const_1 

        #try:
        sol = fsolve(func, x0 = np.array([d1_0, d2_0, d3_0]), args = d)
        if type(sol) == np.ndarray:
            sol = sol.tolist()
        print 'solution: ', sol
        b4 = _sage_const_2  - (sol[_sage_const_0 ] + sol[_sage_const_1 ] + sol[_sage_const_2 ]) 
        if sol[_sage_const_0 ] > _sage_const_1  or sol[_sage_const_1 ] > _sage_const_1  or sol[_sage_const_2 ] > _sage_const_1  or b4 > _sage_const_1  or sol[_sage_const_0 ] < _sage_const_0  or sol[_sage_const_1 ] < _sage_const_0  or sol[_sage_const_2 ] < _sage_const_0  or b4 < _sage_const_0 :
            print 'solution ', [sol, b4], ' is out of bounds.'
            print 'wrong vol: ', Series(b1 = sol[_sage_const_0 ], b2 = sol[_sage_const_1 ], b3 = sol[_sage_const_2 ])
            continue
        #except:
        #    continue
        
        const, vol = constraint(Series, sol, SIDE_LENGTH, out_file)
        

    print 'Done.'

    return vol, sol

def grid_NSolve(Series, d, num_sec, SIDE_LENGTH, out_file):
    length = _sage_const_1 /float(num_sec)
    counter = _sage_const_0 
    vol_list = []
    sol_list = []
    const = _sage_const_0 
    count = _sage_const_0 
    
    for i in range(num_sec):
        for j in range(num_sec):
            for k in range(num_sec):
                print 'grid number: ', counter+_sage_const_1 
                out_file.write('grid #%d\n' % (counter+_sage_const_1 ))
                counter += _sage_const_1 
                b1_range = [i*length,(i+_sage_const_1 )*length]
                b2_range = [j*length,(j+_sage_const_1 )*length]
                b3_range = [k*length,(k+_sage_const_1 )*length]
                bound = [b1_range, b2_range, b3_range]
                print 'bounds: '
                print bound
                out_file.write('grid bound: %s\n' % (str(bound)))
                vol,sol = NSolve(Series, d, SIDE_LENGTH, bound, out_file)
                try:
                    if sol == -_sage_const_1 :
                        print 'grid ', counter, ' dose not work.'
                        continue
                except:
                    if vol != -_sage_const_1 :
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
    for i in range(_sage_const_3 , max_sec+_sage_const_1 ):
        print 'Try grid with ', i, ' length segments.'
        out_file.write('grid with length #%d\n' % i)
        vol_list, sol_list = grid_NSolve(Series, d, i, SIDE_LENGTH, out_file)
        if len(vol_list) != _sage_const_0  and len(sol_list) != _sage_const_0 :
            print 'vol_list ', i, ': ', vol_list
            vol_min = min(vol_list)
            sol_min = sol_list[np.argmin(vol_list)]
            if type(sol_min) == np.ndarray:
                sol_min = sol_min.tolist()
            print 'vol min: ', vol_min
            vol_min_list.append(vol_min)
            sol_min_list.append(sol_min)
    if len(vol_min_list) != _sage_const_0  and len(sol_min_list) != _sage_const_0 :
        print 'vol_min_list: ', vol_min_list
        print 'sol_min_list: ', sol_min_list
        vol_min = min(vol_min_list)
        sol_min = sol_min_list[np.argmin(vol_min_list)]
        return vol_min, sol_min
    else:
        return -_sage_const_1 , -_sage_const_1 

SIDE_LENGTH = _sage_const_5 
b1 = var('b1')
b2 = var('b2')
b3 = var('b3')
b4 = var('b4')
'''
series = 1/((8*b1 + 5*b2 + 6*b3 - 10)*(7*b1 + 6*b2 + 6*b3 - 10)*(b1 + 2*b2 + b3 - 2)*b1) + 1/((6*b1 + 4*b2 + 5*b3 - 8)*(5*b1 + 5*b2 + 5*b3 - 8)*(b1 + 2*b2 + b3 - 2)*b1) + 1/((4*b1 + 3*b2 + 4*b3 - 6)*(3*b1 + 4*b2 + 4*b3 - 6)*(b1 + 2*b2 + b3 - 2)*b1) + 1/((2*b1 + 2*b2 + 3*b3 - 4)*(b1 + 3*b2 + 3*b3 - 4)*(b1 + 2*b2 + b3 - 2)*b1) + 1/((7*b1 + 6*b2 + 6*b3 - 10)*(6*b1 + 4*b2 + 5*b3 - 8)*(b1 - b2)*b1) + 1/((5*b1 + 5*b2 + 5*b3 - 8)*(4*b1 + 3*b2 + 4*b3 - 6)*(b1 - b2)*b1) + 1/((3*b1 + 4*b2 + 4*b3 - 6)*(2*b1 + 2*b2 + 3*b3 - 4)*(b1 - b2)*b1) + 1/((2*b1 - b3)*(b1 - b2)*(b1 - 2*b2 - 2*b3 + 2)*b1) + 1/((b1 + 3*b2 + 3*b3 - 4)*(b1 - b2)*b1*(b2 + 2*b3 - 2)) - 1/((b1 + 2*b2 + b3 - 2)*(b1 - 2*b2 - 2*b3 + 2)*b1*(b2 + 2*b3 - 2)) - 1/((10*b1 + 3*b2 + 6*b3 - 10)*(8*b1 + 2*b2 + 5*b3 - 8)*(b1 - b2)*b2) - 1/((8*b1 + 2*b2 + 5*b3 - 8)*(6*b1 + b2 + 4*b3 - 6)*(b1 - b2)*b2) - 1/((6*b1 + b2 + 4*b3 - 6)*(4*b1 + 3*b3 - 4)*(b1 - b2)*b2) - 1/((4*b1 + 3*b3 - 4)*(2*b1 - b2 + 2*b3 - 2)*(b1 - b2)*b2) + 1/((2*b1 - b2 + 2*b3 - 2)*(b1 - b2)*(2*b2 - b3)*b2) - 1/((4*b1 + 4*b2 + b3 - 4)*(2*b1 - b3)*(2*b2 - b3)*b3)
'''
series = -_sage_const_1 /_sage_const_2 /((_sage_const_4 *b1 + _sage_const_4 *b2 + _sage_const_3 *b3 - _sage_const_4 )*(b1 + b2 + b3 - _sage_const_1 )*b1*b2) - _sage_const_1 /_sage_const_2 /((b1 + b2 + b3 - _sage_const_1 )*b1*b2*b3)

h1 = _sage_const_0 
h2 = _sage_const_2 
h3 = _sage_const_0 
prism, series = Triang_prism(h1,h2,h3)
print 'series for: ', h1, h2, h3
print series
out_path = '%d_%d_%d_test.txt' % (h1,h2,h3)
out_file = open(out_path, 'w')
vol, sol = iter_grid_NSolve(series, _sage_const_4 , SIDE_LENGTH, out_file)
out_file.close()

