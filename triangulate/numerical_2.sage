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
    
    #print 'Hilb: ', str(Hilb(t1 = (m*b1).exp(), t2 = (m*b2).exp(), t3 = (m*b3).exp(), t4 = (m*4).exp())).replace('e', 'E')
    
    Series = Hilb(t1 = (m*b1).exp(), t2 = (m*b2).exp(), t3 = (m*b3).exp(), t4 = (m*4).exp()).series(m==0, 1)
    
    Series = Series.truncate()
    
    return Series

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

# Triangulate the square prism
def init_square_stack(h):
    if h == 0:
        prism = []
        series = 0
        return prism, series
    # Initalize the cube
    cube_1 = [[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]]
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
    if h == 1:
        series = Hilb(triang)
        return triang, series
    for z in range(1, h):
        move = np.array([0,0,z])
        for i in range(6):
            # The i-th tetrahedron
            tetra_moved = []
            for j in range(4):
                # The j-th vertex
                vert_moved = (np.array(triang_cube_1[i][j])+move).tolist()
                tetra_moved.append(vert_moved)
            triang.append(tetra_moved)
    
    hilb = Hilb(triang)
    
    return triang, hilb

def Triang_square(h1,h2,h3,h4, orient):
    
    # Input: ordered heights (1,0,h1), (0,0,h2), (0,1,h2), (1,1,h4)
    # Find the Orientations:
    # Two orientations (trans:0 or cis:1)
    h_list = [h1, h2, h3, h4]
    h_min = min(h_list)
    h_list.remove(h_min)
    h_min_2 = min(h_list)
    h_list.remove(h_min_2)
    h_min_3 = min(h_list)
    h_list.remove(h_min_3)
    h_max = h_list[0]
    assert h_max >= h_min_3 >= h_min_2 >= h_min
    
    # Phase 1: Cube blocks up to h_min
    cube, hilb = init_square_stack(h_min)
    # Phase 2: Tetra transition to triangle base
    if orient == 0:
        # 1243
        for i in range(h_min_2 - h_min):
            prism_1 = [[[1,0,h_min+i],[1,1,h_min],[0,1,h_min+i],[1,0,h_min+i+1]]]
            prism_2 = [[[1,0,h_min+i+1],[1,1,h_min],[0,1,h_min+i],[0,1,h_min+i+1]]]
            prism_3 = [[[0,0,h_min+i],[1,0,h_min+i],[0,1,h_min+i],[0,0,h_min+i+1]]]
            prism_4 = [[[1,0,h_min+i],[1,0,h_min+i+1],[0,0,h_min+i+1],[0,1,h_min+i]]]
            prism_5 = [[[1,0,h_min+i+1],[0,1,h_min+i],[0,1,h_min+i+1],[0,0,h_min+i+1]]]
            prism = prism_1+prism_2+prism_3+prism_4+prism_5
            cube += prism
            
            series_1 = Hilb(prism_1)
            series_2 = Hilb(prism_2)
            series_3 = Hilb(prism_3)
            series_4 = Hilb(prism_4)
            series_5 = Hilb(prism_5)
            hilb += series_1 + series_2 + series_3 + series_4 + series_5
        
        for i in range(h_min_3 - h_min_2):
            prism_1 = [[[1,0,h_min_2],[1,1,h_min],[0,0,h_min_2+i],[0,0,h_min_2+i+1]]]
            prism_2 = [[[0,0,h_min_2+i],[1,1,h_min],[0,1,h_min_2+i],[0,0,h_min_2+i+1]]]
            prism_3 = [[[0,0,h_min_2+i+1],[0,1,h_min_2+i+1],[0,1,h_min_2+i],[1,1,h_min]]]
            prism = prism_1+prism_2+prism_3
            cube += prism
            
            series_1 = Hilb(prism_1)
            series_2 = Hilb(prism_2)
            series_3 = Hilb(prism_3)
            hilb += series_1 + series_2 + series_3
        
        for i in range(h_max - h_min_3):
            prism_1 = [[[1,0,h_min_2], [1,1,h_min],[0,0,h_min_3+i],[0,0,h_min_3+i+1]]]
            prism_2 = [[[0,1,h_min_3],[1,1,h_min],[0,0,h_min_3+i],[0,0,h_min_3+i+1]]]
            prism = prism_1+prism_2
            cube += prism
            
            series_1 = Hilb(prism_1)
            series_2 = Hilb(prism_2)
            hilb += series_1 + series_2
        
    if orient == 1:
        # 1234
        for i in range(h_min_2 - h_min):
            prism_1 = [[[1,0,h_min],[1,1,h_min+i],[1,1,h_min+i+1],[0,0,h_min+i+1]]]
            prism_2 = [[[1,0,h_min],[1,1,h_min+i],[0,0,h_min+i],[0,0,h_min+i+1]]]
            prism_3 = [[[0,0,h_min+i],[1,1,h_min+i],[0,1,h_min+i],[0,0,h_min+i+1]]]
            prism_4 = [[[0,0,h_min+i+1],[1,1,h_min+i],[0,1,h_min+i],[1,1,h_min+i+1]]]
            prism_5 = [[[0,0,h_min+i+1],[1,1,h_min+i+1],[0,1,h_min+i],[0,1,h_min+i+1]]]
            prism = prism_1+prism_2+prism_3+prism_4+prism_5
            cube += prism
            
            series_1 = Hilb(prism_1)
            series_2 = Hilb(prism_2)
            series_3 = Hilb(prism_3)
            series_4 = Hilb(prism_4)
            series_5 = Hilb(prism_5)
            hilb += series_1 + series_2 + series_3 + series_4 + series_5
        
        for i in range(h_min_3 - h_min_2):
            prism_1 = [[[1,0,h_min],[1,1,h_min],[0,0,h_min_2+i],[0,0,h_min_2+i+1]]]
            prism_2 = [[[0,0,h_min_2+i],[1,1,h_min],[0,1,h_min_2+i],[0,0,h_min_2+i+1]]]
            prism_3 = [[[0,0,h_min_2+i+1],[0,1,h_min_2+i+1],[0,1,h_min_2+i],[1,1,h_min]]]
            prism = prism_1+prism_2+prism_3
            cube += prism
            
            series_1 = Hilb(prism_1)
            series_2 = Hilb(prism_2)
            series_3 = Hilb(prism_3)
            hilb += series_1 + series_2 + series_3
        
        for i in range(h_max - h_min_3):
            prism_1 = [[[1,0,h_min], [1,1,h_min],[0,0,h_min_3+i],[0,0,h_min_3+i+1]]]
            prism_2 = [[[0,1,h_min_3],[1,1,h_min],[0,0,h_min_3+i],[0,0,h_min_3+i+1]]]
            prism = prism_1+prism_2
            cube += prism
            
            series_1 = Hilb(prism_1)
            series_2 = Hilb(prism_2)
            hilb += series_1 + series_2
    if orient == 2:
        # 1423
        for i in range(h_min_2 - h_min):
            prism_1 = [[[1,0,h_min+i],[1,1,h_min+i],[1,1,h_min+i+1],[0,0,h_min+i+1]]]
            prism_2 = [[[1,0,h_min+i],[1,1,h_min+i],[0,0,h_min+i],[0,0,h_min+i+1]]]
            prism_5 = [[[0,0,h_min+i+1],[1,1,h_min+i+1],[1,0,h_min+i],[0,1,h_min]]]
            prism_3 = [[[0,0,h_min+i],[1,1,h_min+i],[0,1,h_min],[0,0,h_min+i+1]]]
            prism_4 = [[[0,0,h_min+i+1],[1,1,h_min+i],[0,1,h_min],[1,1,h_min+i+1]]]
            
            prism = prism_1+prism_2+prism_3+prism_4+prism_5
            cube += prism
            
            series_1 = Hilb(prism_1)
            series_2 = Hilb(prism_2)
            series_3 = Hilb(prism_3)
            series_4 = Hilb(prism_4)
            series_5 = Hilb(prism_5)
            hilb += series_1 + series_2 + series_3 + series_4 + series_5
        
        
        for i in range(h_min_3 - h_min_2):
            
            prism_1 = [[[1,0,h_min_2],[0,0,h_min_2+i],[1,1,h_min_2+i],[0,0,1+h_min_2+i]]]
            prism_2 = [[[1,0,h_min_2],[1,1,h_min_2+i],[1,1,1+h_min_2+i],[0,0,1+h_min_2+i]]]
            prism_3 = [[[0,0,h_min_2+i],[1,1,h_min_2+i],[0,1,h_min],[0,0,1+h_min_2+i]]]
            prism_4 = [[[1,1,1+h_min_2+i],[1,1,h_min_2+i],[0,1,h_min],[0,0,1+h_min_2+i]]]
            prism = prism_1+prism_2+prism_3+prism_4
            cube += prism
            
            series_1 = Hilb(prism_1)
            series_2 = Hilb(prism_2)
            series_3 = Hilb(prism_3)
            series_4 = Hilb(prism_4)
            hilb += series_1 + series_2 + series_3 + series_4
        
        
        for i in range(h_max - h_min_3):
            prism_1 = [[[1,0,h_min_2],[1,1,h_min_3],[0,0,h_min_3+i],[0,0,1+h_min_3+i]]]
            prism_2 = [[[0,1,h_min],[1,1,h_min_3],[0,0,h_min_3+i],[0,0,1+h_min_3+i]]]
            prism = prism_1+prism_2
            cube += prism
            
            series_1 = Hilb(prism_1)
            series_2 = Hilb(prism_2)
            hilb += series_1 + series_2
            
    return cube, hilb

def generate_triang_prism(max_height, num_height):
    for height in range(1, max_height+1):
        print 'Height: ', height
        SIDE_LENGTH = int((height+1)/2)
        if num_height <= 1:
            print('Wrong input')
            return -1
        if num_height == 2:
            for h1 in range(1, height+1):
                h2 = height-h1
                if h2 > 0 and h2 <= h1:
                    prism, series = Triang_prism(0, h1, h2)
                    print 'param: ', [0, h1, h2]
                    print 'prism: ', prism
                    print 'series: ', series
        if num_height == 3:
            for h1 in range(1, N+1):
                for h2 in range(1, h1+1):
                    h3 = N-h1-h2
                    if h3 > 0 and h3 <= h2:
                        prism, series = Triang_prism(h1, h2, h3)
                        print 'param: ', [h1, h2, h3]
                        print 'prism: ', prism
                        print 'series: ', series
    print 'Done.'

def check_tri_prism(input_path, output_path):
    input_file = open(input_path, 'r')
    output_file = open(output_path, 'w')
    vol_dict = {}
    counter = 0
    max_count = 100
    for line in input_file:
        if counter >= max_count:
            break
        data = eval(line)
        x = float(data[0][0])
        y = float(data[0][1])
        z = float(data[0][2])
        vol = float(data[1])
        if vol <= 0:
            continue
        key = '%f_%f_%f' % (x,y,z)
        try:
            val = vol_dict[key]
            if abs(val - vol) <= 0.000001:
                continue
            else:
                output_file.write('vol:%f\n'%(vol))
                print 'param: ', x,y,z
                print 'different vol: ', vol
        except:
            vol_dict[key] = vol
            prism, series = Triang_prism(int(x), int(y), int(z))
            output_file.write('\nparam: [%f,%f,%f]\nprism: %s\nseries: %s\nvol: %f\n' % (x,y,z, str(prism), str(series), vol))
        counter += 1
    input_file.close()
    output_file.close()
    print 'Done.'
        
def check_sq_prism(input_path, output_path):
    input_file = open(input_path, 'r')
    output_file = open(output_path, 'w')
    vol_dict = {}
    counter = 0
    max_count = 100
    orient = 0
    for line in input_file:
        if counter >= max_count:
            break
        data = eval(line)
        x = float(data[0][0])
        y = float(data[0][1])
        z = float(data[0][2])
        w = float(data[0][3])
        vol = float(data[1])
        if vol <= 0:
            continue
        key = '%f_%f_%f_%f' % (x,y,z,w)
        try:
            val = vol_dict[key]
            if abs(val - vol) <= 0.000001:
                continue
            else:
                orient = (orient + 1)%2
                output_file.write('orient: %d, vol:%f\n'%(orient, vol))
                print 'param: ', x,y,z,w, orient
                print 'different vol: ', vol
        except:
            vol_dict[key] = vol
            prism, series = Triang_square(int(x), int(y), int(z), int(w), 0)
            output_file.write('\nparam: [%d,%d,%d,%d]\nprism: %s\nseries: %s\norient: %d, vol: %f\n' 
                              % (int(x),int(y),int(z),int(w), str(prism), str(series), 0, vol))
        counter += 1
    input_file.close()
    output_file.close()
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

def func(p, *d):
    f1, f2, f3 = d
    return (f1(b1 = p[0], b2 = p[1], b3 = p[2]), f2(b1 = p[0], b2 = p[1], b3 = p[2]), f3(b1 = p[0], b2 = p[1], b3 = p[2]))

def test_func_1(x, a):
    return a * float(1.0)/(x)
def test_func_2(x, a, b, c):
    return -1*a * np.log(x-b) + c
def test_func(x):
    return test_func_1((x+10), 4.62386348)+0.03
def test_func_inv(y):
    return 4.62386348/(y-0.03)-10
'''
def dist_to_func(pt):
    x0 = float(pt[0])
    y0 = 500*float(pt[1])
    bnds = np.array([(0,None)])
    dist =  lambda x: np.sqrt((x-x0)**2 + (500*test_func(x)-y0)**2)
    res = minimize(dist, [1], bounds=bnds)
    print res
    return res

def dist_to_func(pt):
    x0 = float(pt[0])
    y0 = float(pt[1])
    if x0 <= 100:
        return abs(x0 - test_func_inv(y0))
    else:
        return abs(y0 - test_func(x0))
'''
def dist_to_func(pt):
    x0 = float(pt[0])
    y0 = float(pt[1])
    return abs(y0 - test_func(x0))
    
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
'''
def constraint(Series, sol, vol_range):
    vol = Series(b1 = sol[0], b2 = sol[1], b3 = sol[2])
    vol = abs(vol)
    vol_min = vol_range[0]
    vol_max = vol_range[1]
    if vol <= 1 and vol >= float(1/((3*SIDE_LENGTH)**3)):
        return 1, vol
    
    print 'volume: ', vol, ' is out of bounds.'
    
    return 0, -1
'''

def constraint(Series, sol):
    global vol_min_global
    vol = Series(b1 = sol[0], b2 = sol[1], b3 = sol[2])
    vol = abs(vol)
    if vol <= 1 and vol >= vol_min_global:
        return 1, vol
    
    print 'volume: ', vol, ' is out of bounds.'
    
    return 0, -1
    
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
        print 'reset starting point: ', d1_0, d2_0, d3_0

        try:
            sol = fsolve(func, x0 = np.array([d1_0, d2_0, d3_0]), args = d)
            print 'solution: ', sol
            print 'guessed vol: ', Series(b1 = sol[0], b2 = sol[1], b3 = sol[2])
        except:
            continue
        
        const, vol = constraint(Series, sol)

    print 'Done.'

    return vol, sol

def grid_NSolve(Series, d, num_sec, vol_range, out_file):
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
                #out_file.write('grid #%d\n' % (counter+1))
                counter += 1
                b1_range = [i*length,(i+1)*length]
                b2_range = [j*length,(j+1)*length]
                b3_range = [k*length,(k+1)*length]
                bound = [b1_range, b2_range, b3_range]
                print 'bounds: '
                print bound
                #out_file.write('grid bound: %s\n' % (str(bound)))
                vol,sol = NSolve(Series, d, vol_range, bound, out_file)
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
    
    return vol_list, sol_list

def fit_NSolve(Series, max_range_start, heights, vol_range):
    d1 = diff(Series, b1)
    d2 = diff(Series, b2)
    d3 = diff(Series, b3)
    d = (d1, d2, d3)
    MAX_NUM_VOL = 10
    target_vol = vol_range[0]
    max_diff = vol_range[1]
    range_dict = {}
    height = max(heights)
    idx = get_vol_idx(heights)
    print 'vol idx: ', idx
    
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
                    print 'try bound: ', bound
                    range_dict['%d_%d_%d' % (i,j,k)] = bound
                    vol, sol = NSolve(Series, d, bound)
                    if type(sol) == int or type(vol) == int:
                        print 'range ', i,j,k,' does not work'
                        continue
                    if type(sol) == np.ndarray:
                        sol = sol.tolist()
                    target_dist = target_vol - vol
                    fit_dist = dist_to_func([idx, vol])
                    print 'dist to fit function: ', fit_dist
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
        print 'Done.'
        vol_list = list(vol_list)
        sol_list = list(sol_list)
        print 'vol_list: ', vol_list
        print 'sol_list: ', sol_list
        print 'dist_list: ', dist_list
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
                    print 'distance ', dist, ' is too far from target'
                    continue
        try:
            if type(sol_ret) != int and len(sol_ret) != 0:
                print 'good'
                return vol_ret, sol_ret
        except:
            continue
    print 'no valid solution'
    return vol_ret, sol_ret

def expand_NSolve(Series, max_range_start, height, vol_range):
    d1 = diff(Series, b1)
    d2 = diff(Series, b2)
    d3 = diff(Series, b3)
    d = (d1, d2, d3)
    MAX_NUM_VOL = 10
    target_vol = vol_range[0]
    max_diff = vol_range[1]
    range_dict = {}
    dist_min_pos = 99999
    dist_min_neg = 99999
    dist_min_out = 99999
    vol_list_pos = []
    sol_list_pos = []
    vol_list_neg = []
    sol_list_neg = []
    vol_list_out = []
    sol_list_out = []
    
    for max_range in range(max_range_start, max_range_start + height):
        vol_list = []
        sol_list = []
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
                    print 'try bound: ', bound
                    range_dict['%d_%d_%d' % (i,j,k)] = bound
                    vol, sol = NSolve(Series, d, bound)
                    try:
                        dist = target_vol - vol
                        if abs(dist) < 1e-6:
                            if type(sol) == np.ndarray:
                                sol = sol.tolist()
                            return vol, sol
                        if sol == -1:
                            print 'range ', i,j,k,' does not work'
                        continue
                    except:
                        if vol != -1:
                            print 'result vol: ', vol
                            print 'result sol: ', sol

                            if type(sol) == np.ndarray:
                                sol = sol.tolist()
                            vol_list.append(vol)
                            sol_list.append(sol)
        print 'Done.'
        print 'vol_list: ', vol_list
        print 'sol_list: ', sol_list
        vol = -1
        sol = -1
        
        if len(vol_list) != 0 and len(sol_list) != 0:
            for i in range(len(vol_list)):
                vol_tmp = vol_list[i]
                sol_tmp = sol_list[i]
                dist = target_vol - vol_tmp
                if abs(dist) < 1e-6:
                    return vol_tmp, sol_tmp
                if abs(dist) > max_diff:
                    if abs(dist) < dist_min_out:
                        print 'distance ', dist, ' is too far from target'
                        dist_min_out = abs(dist)
                        vol_list_out = vol_tmp
                        sol_list_out = sol_tmp
                    continue
                if dist > 0:
                    if abs(dist) < dist_min_pos:
                        print 'vol is smaller than target'
                        dist_min_pos = abs(dist)
                        vol_list_pos = vol_tmp
                        sol_list_pos = sol_tmp
                    continue
                if dist < 0:
                    if abs(dist) < dist_min_neg:
                        print 'vol is greater than target'
                        dist_min_neg = abs(dist)
                        vol_list_neg = vol_tmp
                        sol_list_neg = sol_tmp
                    continue
        print 'vol_list_pos: ', vol_list_pos
        print 'sol_list_pos: ', sol_list_pos
        print 'vol_list_neg: ', vol_list_neg
        print 'sol_list_neg: ', sol_list_neg
        print 'vol_list_out: ', vol_list_out
        print 'sol_list_out: ', sol_list_out
        print 'max_range: ', max_range
        if len(sol_list_pos) != 0:
            print 'good'
            return vol_list_pos, sol_list_pos
        if len(sol_list_neg) != 0:
            print 'meh'
            if max_range >= max_range_start + height - 1:
                return vol_list_neg, vol_list_neg
            continue
        if len(sol_list_out) != 0:
            print 'nah'
            if max_range >= max_range_start + height - 1:
                return vol_list_out, vol_list_out
            continue
    return -1, -1

def generate_series(min_height, max_height, out_path):
    global vol_min_global
    vol_dict = {}
    vol_list = []
    num_vol = 0
    for h1 in range(min_height, max_height+1):
        target_vol = 16.0/27/h1
        vol_min_global = 1.0/(h1+1)**3
        for h2 in range(0, h1+1):
            h2 = h1 - h2
            for h3 in range(0, h2+1):
                h3 = h2 - h3
                try:
                    data_tmp = vol_dict['%d_%d_%d' % (h1,h2,h3)]
                    continue
                except:
                    pass
                print h1,h2,h3
                prism, series = lift_prism(h1,h2,h3)
                print 'series: ', series
                out_file = open(out_path, 'a')
                out_file.write('[%s,%s]\n' % (str([h1,h2,h3]), str(series)))
                out_file.close()
                num_vol += 1
                if num_vol > 10:
                    print 'Done.'
                    return
                
def generate_vol_2(min_height, max_height, out_path, max_range=0.02):
    global vol_min_global
    vol_dict = {}
    vol_list = []
    num_vol = 0
    for h1 in range(min_height, max_height+1):
        target_vol = 16.0/27/h1
        vol_min_global = 1.0/(h1+1)**3
        for h2 in range(0, h1+1):
            h2 = h1 - h2
            for h3 in range(0, h2+1):
                h3 = h2 - h3
                try:
                    data_tmp = vol_dict['%d_%d_%d' % (h1,h2,h3)]
                    continue
                except:
                    pass
                print h1,h2,h3
                prism, series = lift_prism(h1,h2,h3)
                vol, sol = fit_NSolve(series, 3, [h1,h2,h3], [target_vol, max_range])
                vol_dict['%d_%d_%d' % (h1,h2,h3)] = [vol, sol]
                vol_list.append(vol)
                out_file = open(out_path, 'a')
                out_file.write('[%s,%s,%s]\n' % (str([h1,h2,h3]), str(vol), str(sol)))
                out_file.close()
                num_vol += 1
                print 'vol: ', vol

def generate_vol(min_height, max_height, max_range, out_path):
    global vol_min_global
    vol_dict = {}
    vol_list = []
    num_vol = 0
    for h1 in range(min_height, max_height+1):
        print h1, 0, 0
        num_vol_prev = (h1+1)*h1/2
        print 'num_vol_prev: ', num_vol_prev
        
        vol_min_global = 1.0/(h1+1)**3

        # we assume that vol(h1,0,0) is going to be
        # smaller than vol(h1-1,0,0), and that it is
        # probably close to 16/27/h1
        try:
            max_vol = vol_list[num_vol_prev]
            print 'max_vol: ', max_vol
        except:
            max_vol = 16.0/27/h1
        min_vol_tmp = 16.0/27/(h1+1.0)
        min_vol = min_vol_tmp
        target_vol = max_vol
        max_range = abs(max_vol - min_vol)
        
        print 'target: ', target_vol, ', range: ', max_range
        prism, series = lift_prism(h1,0,0)
        vol, sol = expand_NSolve(series, 5, h1, [target_vol, max_range])
        print 'vol: ', vol
        vol_dict['%d_%d_%d' % (h1,0,0)] = [vol, sol]
        vol_list.append(vol)
        num_vol += 1
        target_vol = vol
        max_range = abs(target_vol - min_vol)

        for h2 in range(0, h1+1):
            h2 = h1 - h2
            for h3 in range(0, h2+1):
                h3 = h2 - h3
                try:
                    data_tmp = vol_dict['%d_%d_%d' % (h1,h2,h3)]
                    continue
                except:
                    pass
                print h1,h2,h3
                print 'target: ', target_vol, ', range: ', max_range
                prism, series = lift_prism(h1,h2,h3)
                vol, sol = expand_NSolve(series, 5, h1, [target_vol, max_range], out_file)
                vol_dict['%d_%d_%d' % (h1,h2,h3)] = [vol, sol]
                vol_list.append(vol)
                num_vol += 1
                print 'vol: ', vol
    print vol_list

    for key, data in vol_dict.items():
        out_file.write(str(key) + ' >>> ' + str(data) + '\n')
    return vol_list

def iter_grid_NSolve(Series, max_sec, SIDE_LENGTH, out_file):
    d1 = diff(Series, b1)
    d2 = diff(Series, b2)
    d3 = diff(Series, b3)
    d = (d1, d2, d3)
    vol_min_list = []
    sol_min_list = []
    for i in range(4, max_sec+1):
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

def clean_err_vol(err_path, out_path):
    global vol_min_global
    err_file = open(err_path, 'r')
    for line in err_file:
        data = eval(line)
        h = data[0]
        vol_err = data[1]
        vol_min_global = 1/(max(h)**3)
        prism, series = lift_prism(h[0],h[1],h[2])
        vol, sol = fit_NSolve(series, 3, h, [16.0/27/max(h), 0.02])
        print vol
        out_file = open(out_path, 'a')
        out_file.write("[%s, %f, %s]\n" % (str(h), vol, sol))
        out_file.close()
    err_file.close()

'''
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

h1 = 2
h2 = 1
h3 = 0
max_height = max([h1,h2,h3])
SIDE_LENGTH = max_height
#SIDE_LENGTH = (max_height+1)/2
prism, series = lift_prism(h1,h2,h3)
print 'series for: ', h1, h2, h3
print series


vol_min_global = 1/(10**3)
#input_min_height = 1
#input_max_height = 1
#print input_min_height, input_max_height
out_path = '12_1_0_test.txt'
generate_vol_2(10, 11, 0.02, out_path)
#out_file = open(out_path, 'w')
#prism, series = lift_prism(10,2,1)
#vol, sol = expand_NSolve(series, 5, 12, [16.0/27/12, abs(16.0/27/12 - 16.0/27/13)], out_file)
#vol, sol = iter_grid_NSolve(series, 5, SIDE_LENGTH, out_file)
#vol_list = generate_vol(input_min_height, input_max_height, 5, out_file)
#vol, sol = expand_NSolve(series, 3, 10, [16.0/27/10, abs(16.0/27/10 - 16.0/27/11)])
#vol, sol = fit_NSolve(series, 3, [10,2,1], [16.0/27/10, 0.02])
#out_file.close()

h_min = 35
h_max = 36
out_path = '/home/ubuntu/Calabi_Yau/triangulate/lift_vol_%d_%d.txt'%(h_min, h_max)
out_file = open(out_path, 'w')
out_file.close()
generate_vol_2(h_min, h_max, out_path=out_path)
print 'completed.'
'''
h_min = 35
h_max = 36
out_path = '/home/ubuntu/Calabi_Yau/triangulate/series_%d_%d.txt'%(h_min, h_max)
out_file = open(out_path, 'w')
out_file.close()
generate_series(h_min, h_max, out_path)