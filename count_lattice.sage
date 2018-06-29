import math
import numpy as np
import scipy
import scipy.optimize
#import matplotlib.pyplot as plt
import csv

def exist(pts, latt):
    latt = np.array(latt)
    for i in range(pts.shape[0]):
        if pts[i][0]==latt[0]:
            if pts[i][1]==latt[1]:
                if pts[i][2]==latt[2]:
                    return 1
    return 0

def dist(p1, p2):
    return sqrt((p1[0]-p2[0])^2+(p1[1]-p2[1])^2+(p1[2]-p2[2])^2)

def on_edge(latt, poly):
    edges = poly.faces(1)
    num_edges = len(edges)
    for i in range(num_edges):
        pt1 = list(edges[i].vertices()[0])
        pt2 = list(edges[i].vertices()[1])
        if (dist(pt1, pt2) == (dist(pt1, latt) + dist(pt2, latt))):
            return 1
    return 0

def on_face(latt, poly):
    faces = poly.faces(2)
    for face in poly.faces(2):
        face_pts = [list(face.vertices()[i]) for i in range(len(face.vertices()))]
        face_poly = Polyhedron(face_pts)
        if face_poly.contains(latt) == 1:
            return 1
    return 0

def count_pts(pts):
    # Count the number of corner points, edge points, face points, and body points
    num_corner = len(pts)
    num_edge = 0
    num_face = 0
    num_body = 0
    #edge = []
    #face = []
    #body = []
    pts_max = int(max(np.amax(pts, axis=0)))+1
    pts_min = int(min(np.amin(pts, axis=0)))-1
    #print 'pts_max: ', pts_max
    #print 'pts_min: ', pts_min
    poly = Polyhedron(pts)
    pts_new = pts
    for i in range(pts_min, pts_max):
        for j in range(pts_min, pts_max):
            for k in range(pts_min, pts_max):
                latt = [i,j,k]
                if exist(np.array(pts), latt) == 1:
                    continue
                if on_edge(latt, poly) == 1:
                    num_edge += 1
                    #edge.append(latt)
                elif on_face(latt, poly) == 1:
                    num_face += 1
                    #face.append(latt)
                elif poly.interior_contains(latt) == 1:
                    num_body += 1
                    #body.append(latt)
    #print 'edge: ', edge
    #print 'face: ', face
    #print 'body: ', body
    return [num_corner, num_edge, num_face, num_body]

input_path = 'output/polygon/cube/'
output_path = 'output/topology/cube/'
file_path = ['1x1.txt', '2x2.txt', '3x3.txt', '4x4.txt', '5x5.txt']
for path in file_path:
    input_file = open(input_path + path , 'r')
    output_file = open(output_path + path, 'w')
    for line in input_file:
        pts = eval(line)
        topo_list = count_pts(pts)
        output_file.write('%s\n' % topo_list)
        print('Done.')


