
# This file was *autogenerated* from the file generate_poly.sage
from sage.all_cmdline import *   # import sage library

_sage_const_3 = Integer(3); _sage_const_2 = Integer(2); _sage_const_1 = Integer(1); _sage_const_0 = Integer(0); _sage_const_2p5 = RealNumber('2.5'); _sage_const_4 = Integer(4); _sage_const_100 = Integer(100); _sage_const_2048 = Integer(2048)
import math
import numpy as np
import scipy
import scipy.optimize
#import matplotlib.pyplot as plt
import csv

#Helper Functions

def exist(pts, latt):
    latt = np.array(latt)
    for i in range(pts.shape[_sage_const_0 ]):
        if pts[i][_sage_const_0 ]==latt[_sage_const_0 ]:
            if pts[i][_sage_const_1 ]==latt[_sage_const_1 ]:
                if pts[i][_sage_const_2 ]==latt[_sage_const_2 ]:
                    return _sage_const_1 
    return _sage_const_0 

def face_list(p1):
    face_p1 = p1.faces(_sage_const_1 )
    list_face_p1 = list(face_p1)
    faces = []
    for i in range(len(face_p1)):
        faces.append([])
    for i in range(len(face_p1)):
        faces[i].append(list(list_face_p1[i].vertices()[_sage_const_0 ]))
        faces[i].append(list(list_face_p1[i].vertices()[_sage_const_1 ]))
    return faces

#Convert a set of vertices into a list
def vert_to_list(vert):
    num_vert = len(vert)
    vert_list = []
    for i in range(num_vert):
        vert_list.append(list(vert[i]))
    return vert_list

#Calcualte the distance between two pts
def dist(p1, p2):
    return sqrt((p1[_sage_const_0 ]-p2[_sage_const_0 ])**_sage_const_2 +(p1[_sage_const_1 ]-p2[_sage_const_1 ])**_sage_const_2 +(p1[_sage_const_2 ]-p2[_sage_const_2 ])**_sage_const_2 )

def on_edge(latt, faces):
    for i in range(len(faces)):
        if (dist(faces[i][_sage_const_0 ], faces[i][_sage_const_1 ]) == (dist(faces[i][_sage_const_0 ], latt) + dist(faces[i][_sage_const_1 ], latt))):
            #print 'edge: ', faces[i][0], ' and ', faces[i][1]
            #print 'l1: ', dist(faces[i][0], faces[i][1])
            #print 'l2: ', dist(faces[i][0], latt) + dist(faces[i][1], latt)
            return _sage_const_1 
    return _sage_const_0 

#Add lattice points onto the edge of a polygon
def add_lattice(poly):
    #convert vertices to points
    pts = []
    vert = list(poly.vertices())
    num_pts = len(vert)

    for i in range(num_pts):
        pts.append(list(vert[i]))

    #find the maximum of points
    pts = np.array(pts)
    pts_max = int(max(np.amax(np.absolute(pts), axis=_sage_const_0 )))+_sage_const_1 
    pts_new = pts

    faces = face_list(poly)

    for i in range(-pts_max, pts_max):
        for j in range(-pts_max, pts_max):
            for k in range(-pts_max, pts_max):
                latt = [i,j,k]
                if latt in pts.tolist():
                    continue
                if poly.contains(latt) == _sage_const_1  or on_edge(latt, faces) == _sage_const_1 :
                #if on_edge(latt, faces) == 1:
                    pts_new = np.append(pts_new, np.array(latt).reshape((_sage_const_1 ,_sage_const_3 )), axis = _sage_const_0 )
    #print 'pts_new: '
    #print pts_new
    pts_new = pts_new.tolist()
    poly_new = Polyhedron(vertices = pts_new)



    return poly_new, pts_new


#remove a point from a polyhedron
#input: points of a polyhedron; a point to be removed
#output: polyhedron with a point removed
def remove_pts(pts, remove_pt):
    #points to remove
    pts_removed = pts

    #backup points
    pts_save = np.array(pts)

    #remove points
    pts_removed.remove(remove_pt)

    #restore
    pts = pts_save.tolist()

    poly_new = Polyhedron(vertices = pts_removed)

    return poly_new

#Count the number of lattice points inside a polytope
def num_latt(poly):
    count = _sage_const_0 

    #convert vertices to points
    pts = []
    vert = list(poly.vertices())
    num_pts = len(vert)

    for i in range(num_pts):
        pts.append(list(vert[i]))

    #find the maximum of points
    pts = np.array(pts)
    pts_max = int(max(np.amax(pts, axis=_sage_const_0 )))+_sage_const_1 

    faces = face_list(poly)

    for i in range(-pts_max, pts_max):
        for j in range(-pts_max, pts_max):
            for k in range(-pts_max, pts_max):
                latt = [i,j,k]
                if latt in pts.tolist():
                    continue
                if poly.contains(latt) == _sage_const_1 :
                    count += _sage_const_1 
    return count




#Execution Cell

output_path = 'output/polygon/poly_out_1024.txt'

def generate_poly(size, num_poly, output_path):

    output = open(output_path, 'w')

    p1 = polytopes.hypercube(_sage_const_3 )
    p1 = p1.translation([_sage_const_1 ,_sage_const_1 ,_sage_const_1 ])
    p1 = p1.dilation(size)

    poly, pts = add_lattice(p1)
    print "Step 1: add lattice done."
    
    for i in range(num_poly):

        face_pts = list(poly.faces(_sage_const_2 ))
        num_faces = len(face_pts)
        face_idx = np.random.randint(num_faces)
        pt_idx = np.random.randint(_sage_const_3 )
        remove_pt = list(face_pts[face_idx].vertices()[pt_idx])
        poly = remove_pts(pts, remove_pt)
        poly_vert = vert_to_list(poly.vertices())
        
        output.write("%s\n" % poly_vert)
        
        #poly.plot().save("img/plot_2_%d.png" % i)
        print 'Polytope ', i, ' done.'

        if len(poly_vert) <= _sage_const_4 :
            print 'Not enough points.'
            break

    print "Step 2: remove points done."
    output.close()

#def tetra():

def cylinder_tri(height):
    out_path = 'output/polygon/cylinder/cyl_tri_%d.txt' % height
    out_file = open(out_path, 'w')
    for N in range(_sage_const_1 , height):
        for h1 in range(_sage_const_1 , N+_sage_const_1 ):
            for h2 in range(_sage_const_1 , h1+_sage_const_1 ):
                h3 = N-h1-h2
                if h3 > _sage_const_0  and h3 <= h2:
                    out = [[_sage_const_0 ,_sage_const_0 ,_sage_const_0 ],[_sage_const_1 ,_sage_const_0 ,_sage_const_0 ],[_sage_const_0 ,_sage_const_1 ,_sage_const_0 ],[_sage_const_0 ,_sage_const_0 ,h1],[_sage_const_1 ,_sage_const_0 ,h2],[_sage_const_0 ,_sage_const_1 ,h3]]
                    out_file.write("%s\n" % out)
    out_file.close()
    print("Done.")

def cylinder_sq(height):
    out_path = 'output/polygon/cylinder/cyl_sq_%d.txt' % height
    out_file = open(out_path, 'w')
    for N in range(_sage_const_1 , height):
        for h1 in range(_sage_const_1 , N+_sage_const_1 ):
            for h2 in range(_sage_const_1 , h1+_sage_const_1 ):
                for h3 in range(_sage_const_1 , h2+_sage_const_1 ):
                    h4 = N-h1-h2-h3
                    if h4 > _sage_const_0  and h4 <= h3:
                        out = [[_sage_const_0 ,_sage_const_0 ,_sage_const_0 ],[_sage_const_1 ,_sage_const_0 ,_sage_const_0 ],[_sage_const_0 ,_sage_const_1 ,_sage_const_0 ],[_sage_const_1 ,_sage_const_1 ,_sage_const_0 ],[_sage_const_0 ,_sage_const_0 ,h1],[_sage_const_1 ,_sage_const_0 ,h2],[_sage_const_0 ,_sage_const_1 ,h3],[_sage_const_1 ,_sage_const_1 ,h4]]
                        out_file.write("%s\n" % out)
    out_file.close()
    print("Done.")

def simplex(height):
    out_path = 'output/polygon/simplex/xyz_%d.txt' % height
    out_file = open(out_path, 'w')
    for x in range(_sage_const_1 , height):
        for y in range(_sage_const_1 , height):
            for z in range(_sage_const_1 , height):
                out = [[x,_sage_const_0 ,_sage_const_0 ],[_sage_const_0 ,y,_sage_const_0 ],[_sage_const_0 ,_sage_const_0 ,z]]
                out_file.write("%s\n" % out)
    out_file.close()
    print("Done.")

'''
size = raw_input("Please input the size of the cube: ")
size = int(size)
num_poly = raw_input("Please input the number of polygon: ")
num_poly = int(num_poly)
output_path = raw_input("Please input the output path: ")
output_path = str(output_path)
'''
size = _sage_const_2p5 
num_poly = _sage_const_2048 
output_path = 'output/polygon/cube/5x5.txt'

#generate_poly(size, num_poly, output_path)
height = _sage_const_100 

cylinder_tri(height)
cylinder_sq(height)
simplex(height)

