import math
import numpy as np
import scipy
import scipy.optimize
from scipy.spatial.distance import euclidean
import matplotlib.image as mpimg
PointConfiguration.set_engine('internal')

def idx_to_pts(triang, pts):
    # Input a list of lists of indicies
    # Output a list of lists of points
    triang_new = []
    for i in range(len(triang)):
        triang_new.append([pts[j] for j in triang[i]])
    return triang_new

def in_list(e, l):
    # Check if a list e is in a list of lists l
    for i in range(len(l)):
        check = 0
        length = len(l[i])
        for j in range(length):
            if e[j] == l[i][j]:
                # if one element is the same, increment
                check += 1
            elif e[j] != l[i][j]:
                break
        if check == length:
            return 1
    return 0

def find_face(p, pts):
    # Given a list of vertices of a polytope, find the faces that contain the point p
    num_pts = len(pts)
    if num_pts == 1:
        P1 = Polyhedron(vertices = [tuple(pts)])
    else:
        P1 = Polyhedron(pts)
    faces = []
    for face in P1.faces(2):
        # face_list is a list of verticies of a face
        face_list = [list(face.vertices()[i]) for i in range(len(face.vertices()))]
        #print 'face_list: ', face_list
        #print 'p: ', p
        if in_list(p, face_list):
            #print 'p in face_list'
            faces.append(face_list)
    #print 'faces: ', faces
    return faces

def find_center(p1, p2, p3):
    #print 'center around: ', p1, p2, p3
    p1 = np.array(p1)
    p2 = np.array(p2)
    p3 = np.array(p3)
    c = (p1+p2+p3)/3
    c = c.tolist()
    return c

def find_adj(pt, face_pts):
    # Input: a reference point, pt, and a list of vertices of a face, face_pts
    # face_pts must have at least 3 elements
    
    # Pick a random point, c, on the plane of the face
    #print 'input face: ', face_pts
    #print 'input point: ', pt
    num_face_pts = len(face_pts)
    if num_face_pts < 3:
        raise ValueError('Face %s has too little points.' % face_pts)
    
    # pick 2 points from the face (arbitrary)
    c = find_center(face_pts[0], face_pts[1], face_pts[num_face_pts - 1])
    c = np.array(c)
    #print 'center: ', c
    pt = np.array(pt)

    # Find the vector, v, from c to the reference point, pt
    v = pt - c
    
    # Find the two points, pt1 and pt2, on that face that is closest to pt
    dist = [euclidean(pt, p) for p in face_pts]
    #print 'dist: ', dist
    len_dist = len(dist)
    stop = 0 # Stop sign to stop the while loop
    
    while stop == 0:
        # Sort the distance list
        dist_sort = sorted(dist)
        #print 'dist_sort: ', dist_sort
        sort_1 = dist_sort[1]
        sort_2 = dist_sort[2]
        if sort_1 == sort_2:
            id1 = dist.index(sort_1)
            dist[id1] = 0
            id2 = dist.index(sort_2)
        else:
            id1 = dist.index(sort_1)
            id2 = dist.index(sort_2)
        pt1 = np.array(face_pts[id1])
        #print 'pt1: ', pt1
        pt2 = np.array(face_pts[id2])
        #print 'pt2: ', pt2

        # Find the vector v1 from c to p1 and v2 from c to p2
        v1 = pt1 - c
        v2 = pt2 - c
        # Calculate cross(v1, v) and cross(v2, v)
        c1 = np.cross(v1, v)
        c2 = np.cross(v2, v)
        n1 = np.linalg.norm(c1)
        n2 = np.linalg.norm(c2)
        c1 = c1/n1
        c2 = c2/n2
        #print 'c1: ', c1
        #print 'c2: ', c2
        # If the two cross products are in opposite directions, then return the 2 points
        diff = np.linalg.norm(c1 + c2)
        if diff <= 1e-10:
            # Different sign
            # If zero, at least one point on the opposite side of pt, still works
            #print 'Opposite'
            stop = 1
            break
        # Else, the two cross products are in the same direction
        # then remove the point that is farther from p from face_pts and dist
        else:
            # Same sign
            #print 'Same'
            #print 'result: ', diff
            dist.remove(sort_2)
            face_pts.remove(pt2.tolist())
            #print 'dist updated: ', dist
            #print 'face_pt updated: ', face_pts
            continue
    pt1 = pt1.tolist()
    pt2 = pt2.tolist()
    return pt1, pt2


def split_tetra(pts):
    # Given the vertices of a polyhedron, split it into two polygon, one of which is a tetrahedron (hopefully)
    for p in pts:
        pts1 = [point for point in pts if point != p]
        # Find all faces that contain p
        #print 'reference point: ', p
        faces = find_face(p, pts)
        num_faces = len(faces)
        if num_faces != 3:
            # If the point is not contained in a tetrahedron, then find another point
            continue
        adj_pts = []
        adj_pts.append(p)
        
        for i in range(num_faces):
            # For each face, find two points that are adjacent to p
            p1, p2 = find_adj(p, faces[i])
            if p1 not in adj_pts:
                adj_pts.append(p1)
            if p2 not in adj_pts:
                adj_pts.append(p2)
        break
    return pts1, adj_pts

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
    #print 'pts_max and pts_min: ', pts_max, pts_min
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

def split_triang(pts):
    # Input points is already checked lattice
    #pts = check_latt(pts)
    print 'input points: ', pts
    #print 'num_pts: ', len(pts)
    pc = PointConfiguration(pts)
    try:
        triang = pc.triangulate()
        triang = idx_to_pts(triang, pts)
    except:
        pts1, pts2 = split_tetra(pts)
        print 'pts1: ', pts1
        print 'pts2: ', pts2
        triang_1 = split_triang(pts1)
        print 'triang1: ', triang_1
        triang_2 = split_triang(pts2)
        print 'triang2: ', triang_2
        triang = triang_1 + triang_2
    return triang

def Triang(pts):
    # Input a list of vertices of a polyhedron
    # Output a complete triangulation containing all lattice points
    pts = check_latt(pts)
    triang = split_triang(pts)
    print 'triangulation: ', triang
    vol_check = 0
    
    while vol_check == 0:
        vol_check = 1
        triang_ret = []
        for tetra in triang:
            print 'tetrahedron: ', tetra
            tet = Polyhedron(tetra)
            vol = tet.volume()
            if vol > 1/6:
                vol_check = 0
                print 'volume ', vol, ' is too large'
                check = 0
                for p in pts:
                    if p in tetra:
                        # If the point is already a corner point
                        print 'corner point, pass.'
                        continue
                    if contain(tet, p) == 1:
                        # If the point is not a corner point but contained inside the tetrahedron
                        print 'contained point ', p, ', append.'
                        check = 1
                        tetra.append(p)
                if check == 1:
                    tetra_triang = split_triang(tetra)
                    print 'updated triang: ', tetra_triang
                    for t in tetra_triang:
                        triang_ret.append(t)
                elif check == 0:
                    triang_ret.append(tetra)
            else:
                triang_ret.append(tetra)
        triang = triang_ret
    
    return triang_ret

out_path = 'output/failed/triang/3x3_test.txt'
out_file = open(out_path, 'w')
pts = [[0.0, 0.0, 3.0], [0.0, 3.0, 0.0], [0.0, 3.0, 3.0], [3.0, 0.0, 0.0], [3.0, 3.0, 0.0], [3.0, 3.0, 3.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 3.0], [3.0, 0.0, 2.0], [3.0, 1.0, 3.0]]
triang = Triang(pts)
out_file.write("%s\n" % triang)
out_file.close()