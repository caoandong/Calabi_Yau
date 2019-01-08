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
    # Add 1 to the end of all vectors
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
    #print 'Hilb: ', Hilb()
    #print Hilb(t1=t, t2=t, t3=t).series(t4, 3)
    #print "p-q web: ", power 
    
    
    m = var('m')
    b1 = var('b1')
    b2 = var('b2')
    b3 = var('b3')
    b4 = var('b4')
    Hilb *= m^4
    
    #print 'Hilb: ', str(Hilb(t1 = (m*b1).exp(), t2 = (m*b2).exp(), t3 = (m*b3).exp(), t4 = (m*(2-b1-b2-b3)).exp())).replace('e', 'E')
    
    
    Series = Hilb(t1 = (m*b1).exp(), t2 = (m*b2).exp(), t3 = (m*b3).exp(), t4 = (m*4).exp()).series(m==0, 1)
    Series = Series.truncate()
    
    return Series

def func(p, *d):
    f1, f2, f3 = d
    return (f1(b1 = p[0], b2 = p[1], b3 = p[2]), f2(b1 = p[0], b2 = p[1], b3 = p[2]), f3(b1 = p[0], b2 = p[1], b3 = p[2]))

def constraint(Series, sol, SIDE_LENGTH):
    vol = Series(b1 = sol[0], b2 = sol[1], b3 = sol[2])
    if vol <= 1 and vol >= float(1/((3*SIDE_LENGTH)**3)):
        return 1, vol

    print 'volume: ', vol, ' is out of bounds.'

    return 0, -1

def NSolve(Series, SIDE_LENGTH):
    d1 = diff(Series, b1)
    d2 = diff(Series, b2)
    d3 = diff(Series, b3)
    d = (d1, d2, d3)
    const = 0
    count = 0

    while const == 0:
        d1_0 = np.random.uniform(low=0, high=5)
        d2_0 = np.random.uniform(low=0, high=5)
        d3_0 = np.random.uniform(low=0, high=5)
        print 'reset starting point: ', d1_0, d2_0, d3_0

        try:
            sol = fsolve(func, x0 = np.array([d1_0, d2_0, d3_0]), args = d)
            print 'solution: ', sol
        except:
            continue
        
        const, vol = constraint(Series, sol, SIDE_LENGTH)

        count += 1
        if count > 1000:
            print 'Infinite loop. Force stop.'
            return -1, -1

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
    if size == 0:
        corner = []
        series = 0
        return corner, series
    # Initalize the cube
    cube_1 = [[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]]
    corner = (size*np.array(cube_1)).tolist()
    # Sample triangulation of a 1x1x1 cube:
    triang_cube_1_1 = [[[0,0,0],[1,0,0],[1,1,0],[0,0,1]]]
    triang_cube_1_2 = [[[0,0,1],[1,0,1],[1,0,0],[1,1,1]]]
    triang_cube_1_3 = [[[0,0,1],[1,1,1],[1,1,0],[1,0,0]]]
    triang_cube_1_4 = [[[0,0,0],[1,1,0],[0,1,0],[0,0,1]]]
    triang_cube_1_5 = [[[0,0,1],[0,1,1],[1,1,0],[1,1,1]]]
    triang_cube_1_6 = [[[0,0,1],[0,1,1],[0,1,0],[1,1,0]]]
    triang_cube_1 = triang_cube_1_1 + triang_cube_1_2 + triang_cube_1_3 + triang_cube_1_4 + triang_cube_1_5 + triang_cube_1_6
    # Stack the cubes together
    triang = triang_cube_1
    if size == 1:
        return triang
    for x in range(0,size):
        for y in range(0,size):
            for z in range(0,size):
                if x==0 and y==0 and z ==0:
                    continue
                move = [x, y, z]
                move = np.array(move)
                for i in range(6):
                    # The i-th tetrahedron
                    tetra_moved = []
                    for j in range(4):
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
    if len(adj_pts) > 3:
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

def Triang_cube(size, num_iteration):
    corner, triang, hilb = init_cube(size)
    # Cut a random corner of the cube
    for i in range(num_iteration):
        idx = np.random.randint(len(corner))
        p = corner[idx]
        print 'i: ', i
        print 'p: ', p
        corner, triang, hilb = cut_corner(p, corner, triang, hilb)
        print 'Num Hilbert: ', len(hilb)
        print 'Num triang: ', len(triang)
        assert len(hilb) == len(triang)
        print ''
        if len(corner) <= 3:
            break

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

def generate_triang_prism(max_height, num_height, train_path, pts_path):
    for height in range(1, max_height+1):
        SIDE_LENGTH = int((height+1)/2)
        if num_height <= 1:
            print('Wrong input')
            return -1
        if num_height == 2:
            for h1 in range(1, height+1):
                h2 = height-h1
                if h2 > 0 and h2 <= h1:
                    prism, series = Triang_prism(0, h1, h2)
                    vol, sol = NSolve(series, SIDE_LENGTH)
                    print 'sol: ', sol
                    print 'vol: ', vol
                    pts = [[0,0,0],[1,0,0],[0,1,0],[0,0,h1],[1,0,h2],[0,1,0]]
                    out = [pts, sol]
                    
                    train_set = []
                    param = [0, h1, h2]
                    print 'param: ', param
                    train_set.append(param)
                    train_set.append(vol)
                    train_file = open(train_path, 'a')
                    pts_file = open(pts_path, 'a')
                    train_file.write("%s\n" % train_set)
                    pts_file.write("%s\n" % out)
                    train_file.close()
                    pts_file.close()
        if num_height == 3:
            for N in range(1, height):
                for h1 in range(1, N+1):
                    for h2 in range(1, h1+1):
                        h3 = N-h1-h2
                        if h3 > 0 and h3 <= h2:
                            prism, series = Triang_prism(h1, h2, h3)
                            vol, sol = NSolve(series, SIDE_LENGTH)
                            
                            
                            print 'sol: ', sol
                            print 'vol: ', vol
                            pts = [[0,0,0],[1,0,0],[0,1,0],[0,0,h1],[1,0,h2],[0,1,h3]]
                            out = [pts, sol]

                            train_set = []
                            param = [h1, h2, h3]
                            print 'param: ', param
                            train_set.append(param)
                            train_set.append(vol)
                            train_file = open(train_path, 'a')
                            pts_file = open(pts_path, 'a')
                            train_file.write("%s\n" % train_set)
                            pts_file.write("%s\n" % out)
                            train_file.close()
                            pts_file.close()
        print 'Height ', height, ' done.'
    print 'Done.'
    
def lift_prism(h1, h2, h3):
    # h1 >= h2 >= h3
    h_list = [h1, h2, h3]
    h_min = min(h_list)
    h_list.remove(h_min)
    h_max = max(h_list)
    h_list.remove(h_max)
    h_mid = h_list[0]
    
    
    if h_min == 0:
        prism = []
        series = 0
        if h_mid > 0:
            for i in range(h_mid):
                prism_1 = [[[1,0,0],[1,1,0],[0,0,i],[0,0,i+1]]]
                prism_2 = [[[0,0,i],[1,1,0],[0,1,i],[0,0,i+1]]]
                prism_3 = [[[1,1,0],[0,1,i],[0,1,i+1],[0,0,i+1]]]
                prism += prism_1 + prism_2 + prism_3

                series_1 = Hilb(prism_1)
                series_2 = Hilb(prism_2)
                series_3 = Hilb(prism_3)
                series += series_1 + series_2 + series_3
        
            for i in range(h_mid, h_max):
                prism_1 = [[[1,0,0],[1,1,0],[0,0,h_mid+i],[0,0,h_mid+i+1]]]
                prism_2 = [[[1,1,0],[0,1,h_mid+i],[0,0,h_mid+i],[0,0,h_mid+i+1]]]
                prism += prism_1 + prism_2

                series_1 = Hilb(prism_1)
                series_2 = Hilb(prism_2)
                series += series_1 + series_2
        elif h_mid == 0:
            for i in range(h_max):
                prism_1 = [[[1,0,0],[1,1,0],[0,0,i],[0,0,i+1]]]
                prism_2 = [[[0,0,i],[1,1,0],[0,1,0],[0,0,i+1]]]
                prism += prism_1 + prism_2
                
                series_1 = Hilb(prism_1)
                series_2 = Hilb(prism_2)
                series += series_1 + series_2
            
        return prism, series
    
    if h_min > 0:
        # Stage 0: base
        prism = [[[1,0,0],[0,0,0],[0,1,0],[1,1,h_min]]]
        series = Hilb(prism)
        # Stage 1: h1 = h2
        for i in range(h_mid):
            prism_1 = [[[1,0,0],[0,0,i],[0,0,i+1],[1,1,h_min]]]
            prism_2 = [[[0,0,i],[0,1,i],[1,1,h_min],[0,0,i+1]]]
            prism_3 = [[[0,1,i],[0,1,i+1],[0,0,i+1],[1,1,h_min]]]
            prism += prism_1 + prism_2 + prism_3
            
            series_1 = Hilb(prism_1)
            series_2 = Hilb(prism_2)
            series_3 = Hilb(prism_3)
            series += series_1 + series_2 + series_3
            
        for i in range(h_mid, h_max):
            prism_1 = [[[1,0,0],[1,1,h_min],[0,0,i],[0,0,i+1]]]
            prism_2 = [[[1,1,h_min],[0,0,i],[0,1,i],[0,0,i+1]]]
            prism += prism_1 + prism_2
            
            series_1 = Hilb(prism_1)
            series_2 = Hilb(prism_2)
            series += series_1 + series_2
    
    return prism, series

def generate_lift_prism(max_height, train_path, pts_path):
    train_file = open(train_path, 'w')
    pts_file = open(pts_path, 'w')
    for height in range(max_height):
        SIDE_LENGTH = int((height+1)/2)
        for h1 in range(1, height+1):
            for h2 in range(0, h1+1):
                h3 = height - h1 - h2
                if h3 >= 0 and h3 <= h2:
                    print 'heights: ', h1, h2 ,h3
                    prism, series = lift_prism(h1,h2,h3)
                    vol, sol = NSolve(series, SIDE_LENGTH)
                    pts = [[0,0,0],[1,0,0],[1,1,h3],[0,1,h2],[0,0,h1]]
                    out = [pts, sol]
                    
                    param = [h1,h2,h3]
                    train_set = [param, vol]
                    
                    train_file.write("%s\n" % train_set)
                    pts_file.write("%s\n" % out)
                    
        print 'height: ', height
    train_file.close()
    pts_file.close()
    print 'Done.'

max_height = 50
#num_height = 2
#SIDE_LENGTH = int((height+1)/2)
train_path = '/home/carnd/CYML/output/train/cylinder/tri_1_to_50_new.txt'
pts_path = '/home/carnd/CYML/output/polygon/cylinder/tri_1_to_50_new.txt'
#train_path = '/home/carnd/CYML/output/train/cylinder/tri_%d_%d.txt' % (height, num_height)
#pts_path = '/home/carnd/CYML/output/polygon/cylinder/tri_%d_%d.txt' % (height, num_height)
#generate_triang_prism(max_height, num_height, train_path, pts_path)
#generate_lift_prism(max_height, train_path, pts_path)
prism, series = Triang_prism(0, 2, 0)

print 'series: '
print series
vol, sol = NSolve(series, 3)
print 'vol: ', vol
print 'sol: ', sol
