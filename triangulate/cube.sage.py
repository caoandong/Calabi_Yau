
# This file was *autogenerated* from the file cube.sage
from sage.all_cmdline import *   # import sage library

_sage_const_3 = Integer(3); _sage_const_2 = Integer(2); _sage_const_1 = Integer(1); _sage_const_0 = Integer(0); _sage_const_6 = Integer(6); _sage_const_5 = Integer(5); _sage_const_4 = Integer(4); _sage_const_100 = Integer(100); _sage_const_1000 = Integer(1000); _sage_const_15 = Integer(15)
import math
import numpy as np
import scipy
import scipy.optimize
import csv
from scipy.optimize import fsolve

PointConfiguration.set_engine('internal')
#PointConfiguration.set_engine('topcom')

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

def dist(p1, p2):
    return sqrt((p1[_sage_const_0 ]-p2[_sage_const_0 ])**_sage_const_2 +(p1[_sage_const_1 ]-p2[_sage_const_1 ])**_sage_const_2 +(p1[_sage_const_2 ]-p2[_sage_const_2 ])**_sage_const_2 )

def on_edge(latt, poly):
    edges = poly.faces(_sage_const_1 )
    num_edges = len(edges)
    for i in range(num_edges):
        pt1 = list(edges[i].vertices()[_sage_const_0 ])
        pt2 = list(edges[i].vertices()[_sage_const_1 ])
        if (dist(pt1, pt2) == (dist(pt1, latt) + dist(pt2, latt))):
            return _sage_const_1 
    return _sage_const_0 

def on_face(latt, poly):
    faces = poly.faces(_sage_const_2 )
    for face in poly.faces(_sage_const_2 ):
        face_pts = [list(face.vertices()[i]) for i in range(len(face.vertices()))]
        face_poly = Polyhedron(face_pts)
        if face_poly.contains(latt) == _sage_const_1 :
            return _sage_const_1 
    return _sage_const_0 

def count_pts(pts):
    # Count the number of corner points, edge points, face points, and body points
    num_corner = len(pts)
    num_edge = _sage_const_0 
    num_face = _sage_const_0 
    num_body = _sage_const_0 
    #edge = []
    #face = []
    #body = []
    pts_max = int(max(np.amax(pts, axis=_sage_const_0 )))+_sage_const_1 
    pts_min = int(min(np.amin(pts, axis=_sage_const_0 )))-_sage_const_1 
    #print 'pts_max: ', pts_max
    #print 'pts_min: ', pts_min
    poly = Polyhedron(pts)
    pts_new = pts
    for i in range(pts_min, pts_max):
        for j in range(pts_min, pts_max):
            for k in range(pts_min, pts_max):
                latt = [i,j,k]
                if exist(np.array(pts), latt) == _sage_const_1 :
                    continue
                if on_edge(latt, poly) == _sage_const_1 :
                    num_edge += _sage_const_1 
                    #edge.append(latt)
                elif on_face(latt, poly) == _sage_const_1 :
                    num_face += _sage_const_1 
                    #face.append(latt)
                elif poly.interior_contains(latt) == _sage_const_1 :
                    num_body += _sage_const_1 
                    #body.append(latt)
    #print 'edge: ', edge
    #print 'face: ', face
    #print 'body: ', body
    return [num_corner, num_edge, num_face, num_body]

print 'Done.'


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
    
    #print 'Hilb: ', str(Hilb(t1 = (m*b1).exp(), t2 = (m*b2).exp(), t3 = (m*b3).exp(), t4 = (m*4).exp())).replace('e', 'E')
    
    
    Series = Hilb(t1 = (m*b1).exp(), t2 = (m*b2).exp(), t3 = (m*b3).exp(), t4 = (m*_sage_const_4 ).exp()).series(m==_sage_const_0 , _sage_const_1 )
    Series = Series.truncate()
    
    return Series



def func(p, *d):
    f1, f2, f3 = d
    return (f1(b1 = p[_sage_const_0 ], b2 = p[_sage_const_1 ], b3 = p[_sage_const_2 ]), f2(b1 = p[_sage_const_0 ], b2 = p[_sage_const_1 ], b3 = p[_sage_const_2 ]), f3(b1 = p[_sage_const_0 ], b2 = p[_sage_const_1 ], b3 = p[_sage_const_2 ]))

def constraint(Series, sol):
    vol = Series(b1 = sol[_sage_const_0 ], b2 = sol[_sage_const_1 ], b3 = sol[_sage_const_2 ])
    if vol <= _sage_const_1  and vol >= float(_sage_const_1 /((_sage_const_3 *SIDE_LENGTH)**_sage_const_3 )):
        return _sage_const_1 , vol

    print 'volume: ', vol, ' is out of bounds.'

    return _sage_const_0 , -_sage_const_1 

def NSolve(Series):
    d1 = diff(Series, b1)
    d2 = diff(Series, b2)
    d3 = diff(Series, b3)
    d = (d1, d2, d3)
    const = _sage_const_0 
    count = _sage_const_0 

    while const == _sage_const_0 :
        d1_0 = np.random.uniform(low=_sage_const_0 , high=_sage_const_5 )
        d2_0 = np.random.uniform(low=_sage_const_0 , high=_sage_const_5 )
        d3_0 = np.random.uniform(low=_sage_const_0 , high=_sage_const_5 )
        print 'reset starting point: ', d1_0, d2_0, d3_0

        try:
            sol = fsolve(func, x0 = np.array([d1_0, d2_0, d3_0]), args = d)
            print 'solution: ', sol
        except:
            continue
        
        const, vol = constraint(Series, sol)

        count += _sage_const_1 
        if count > _sage_const_1000 :
            print 'Infinite loop. Force stop.'
            return -_sage_const_1 , -_sage_const_1 

    print 'Done.'

    return vol, sol


def idx_to_pts(triang, pts):
    # Input a list of lists of indicies
    # Output a list of lists of points
    triang_new = []
    for i in range(len(triang)):
        triang_new.append([pts[j] for j in triang[i]])
    return triang_new

def init_cube(size):
    if size == _sage_const_0 :
        corner = []
        triang = []
        hilb = _sage_const_0 
        return corner, triang, hilb
    # Initalize the cube
    cube_1 = [[_sage_const_0 ,_sage_const_0 ,_sage_const_0 ],[_sage_const_1 ,_sage_const_0 ,_sage_const_0 ],[_sage_const_1 ,_sage_const_1 ,_sage_const_0 ],[_sage_const_0 ,_sage_const_1 ,_sage_const_0 ],[_sage_const_0 ,_sage_const_0 ,_sage_const_1 ],[_sage_const_1 ,_sage_const_0 ,_sage_const_1 ],[_sage_const_1 ,_sage_const_1 ,_sage_const_1 ],[_sage_const_0 ,_sage_const_1 ,_sage_const_1 ]]
    corner = (size*np.array(cube_1)).tolist()
    # Sample triangulation of a 1x1x1 cube:
    triang_cube_1_1 = [[[_sage_const_0 ,_sage_const_0 ,_sage_const_0 ],[_sage_const_1 ,_sage_const_0 ,_sage_const_0 ],[_sage_const_1 ,_sage_const_1 ,_sage_const_0 ],[_sage_const_0 ,_sage_const_0 ,_sage_const_1 ]]]
    triang_cube_1_2 = [[[_sage_const_0 ,_sage_const_0 ,_sage_const_1 ],[_sage_const_1 ,_sage_const_0 ,_sage_const_1 ],[_sage_const_1 ,_sage_const_0 ,_sage_const_0 ],[_sage_const_1 ,_sage_const_1 ,_sage_const_1 ]]]
    triang_cube_1_3 = [[[_sage_const_0 ,_sage_const_0 ,_sage_const_1 ],[_sage_const_1 ,_sage_const_1 ,_sage_const_1 ],[_sage_const_1 ,_sage_const_1 ,_sage_const_0 ],[_sage_const_1 ,_sage_const_0 ,_sage_const_0 ]]]
    triang_cube_1_4 = [[[_sage_const_0 ,_sage_const_0 ,_sage_const_0 ],[_sage_const_1 ,_sage_const_1 ,_sage_const_0 ],[_sage_const_0 ,_sage_const_1 ,_sage_const_0 ],[_sage_const_0 ,_sage_const_0 ,_sage_const_1 ]]]
    triang_cube_1_5 = [[[_sage_const_0 ,_sage_const_0 ,_sage_const_1 ],[_sage_const_0 ,_sage_const_1 ,_sage_const_1 ],[_sage_const_1 ,_sage_const_1 ,_sage_const_0 ],[_sage_const_1 ,_sage_const_1 ,_sage_const_1 ]]]
    triang_cube_1_6 = [[[_sage_const_0 ,_sage_const_0 ,_sage_const_1 ],[_sage_const_0 ,_sage_const_1 ,_sage_const_1 ],[_sage_const_0 ,_sage_const_1 ,_sage_const_0 ],[_sage_const_1 ,_sage_const_1 ,_sage_const_0 ]]]
    triang_cube_1 = triang_cube_1_1 + triang_cube_1_2 + triang_cube_1_3 + triang_cube_1_4 + triang_cube_1_5 + triang_cube_1_6
    # Stack the cubes together
    triang = triang_cube_1
    if size == _sage_const_1 :
        hilb = []
        for tetra in triang:
            series = Hilb([tetra])
            hilb.append(series)
        return corner, triang, hilb
    for x in range(_sage_const_0 ,size):
        for y in range(_sage_const_0 ,size):
            for z in range(_sage_const_0 ,size):
                if x==_sage_const_0  and y==_sage_const_0  and z ==_sage_const_0 :
                    continue
                move = [x, y, z]
                move = np.array(move)
                for i in range(_sage_const_6 ):
                    # The i-th tetrahedron
                    tetra_moved = []
                    for j in range(_sage_const_4 ):
                        # The j-th vertex
                        vert_moved = (np.array(triang_cube_1[i][j])+move).tolist()
                        tetra_moved.append(vert_moved)
                    triang.append(tetra_moved)
    # Also initialize the Hilbert series for each tetrahedron
    hilb = []
    for tetra in triang:
        series = Hilb([tetra])
        hilb.append(series)
    #print 'Init cube'
    #print 'Hilbert series: ', hilb
    return corner, triang, hilb

def cut_corner(p, corner, triang, hilb):
    corner = check_latt(corner)
    corner.remove(p)
    cube_new = Polyhedron(corner)
    vertices = cube_new.vertices()
    corner = [list(vertices[i]) for i in range(len(vertices))]
    #print 'corner: ', corner
    # Find all the tetrahedron points that contain p
    adj_pts = []
    triang_new = []
    hilb_new = []
   
    for tetra in triang:
        if p in tetra:
            # Hilbert
            for vertex in tetra:
                if (vertex != p) and (vertex not in adj_pts):
                    adj_pts.append(vertex)
            
            # I do not trust indexing
            series = Hilb([tetra])
            #print 'remove hilb: ', series
            hilb.remove(series)
        else:
            triang_new.append(tetra)
    triang = triang_new
    #print 'adjacent points: ', adj_pts
    if len(adj_pts) > _sage_const_3 :
        patch = Polyhedron(adj_pts)
        patch_triang = PointConfiguration(adj_pts).triangulate()
        patch_triang = idx_to_pts(patch_triang, adj_pts)
        #print 'patch: ', patch_triang
        for tetra in patch_triang:
            if tetra not in triang:
                triang.append(tetra)
            # Also find the Hilbert Series for each tetrahedron
            series = Hilb([tetra])
            #print 'Patch hilbert: ', series
            hilb.append(series)
    #print 'number of tetrahedron: ', len(triang)
    return corner, triang, hilb

def plot_poly(pts):
    P = Polyhedron(pts)
    P.plot().save("/home/carnd/CYML/img/triang_cube_%d.png" % _sage_const_1 )
    img1 = mpimg.imread("/home/carnd/CYML/img/triang_cube_%d.png" % _sage_const_1 )
    plt.figure(figsize=(_sage_const_15 , _sage_const_15 ))
    plt.imshow(img1)
    plt.show()

def Triang_cube(size, num_iteration):
    corner, triang, hilb = init_cube(size)
    # Cut a random corner of the cube
    series_sum = _sage_const_0 
    for series in hilb:
        series_sum += series
    hilb_ret = [series_sum]
    count_ret = [count_pts(corner)]
    pts_ret = [corner]
    for i in range(num_iteration):
        idx = np.random.randint(len(corner))
        p = corner[idx]
        print 'i: ', i
        print 'p: ', p
        corner, triang, hilb = cut_corner(p, corner, triang, hilb)
        assert len(hilb) == len(triang)
        
        pts_ret.append(corner)
        
        count = count_pts(corner)
        print 'count: ', count
        count_ret.append(count)
        
        series_sum = _sage_const_0 
        for series in hilb:
            series_sum += series
        hilb_ret.append(series_sum)
        
        #plot_poly(corner)
        print ''
        if len(corner) <= _sage_const_3 :
            break
    
    assert len(pts_ret) == len(count_ret) == len(hilb_ret)
    return hilb_ret, pts_ret, count_ret

def vol_cube(size, num_iteration, train_path, count_path, pts_path):
    hilb, pts, count = Triang_cube(size, num_iteration)
    for i in range(len(hilb)):
        vol, sol = NSolve(hilb[i])
        sol = np.around(sol, decimals=_sage_const_4 ).tolist()
        print 'vol: ', vol
        print 'sol: ', sol
        train_set = [pts[i], vol]
        print 'train: ', train_set
        count_set = [count[i], vol]
        print 'count: ', count_set
        pts_set = [pts[i], sol]
        print 'pts: ', pts_set
        
        train_file = open(train_path, 'w')
        count_file = open(count_path, 'w')
        pts_file = open(pts_path, 'w')
        train_file.write("%s\n" % train_set)
        count_file.write("%s\n" % count_set)
        pts_file.write("%s\n" % pts_set)
        train_file.close()
        count_file.close()
        pts_file.close()
    
    print 'Done.'


size = _sage_const_3 
num_iteration = _sage_const_100 
SIDE_LENGTH = size
train_path = '/home/carnd/CYML/output/train/cube/cube_%dx%d_%d.txt' % (size, size, num_iteration)
count_path = '/home/carnd/CYML/output/train/cube/count_%dx%d_%d.txt' % (size, size, num_iteration)
pts_path = '/home/carnd/CYML/output/polygon/cube/pts_%dx%d_%d.txt' % (size, size, num_iteration)
vol_cube(size, num_iteration, train_path, count_path, pts_path)

