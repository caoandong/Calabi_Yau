"""
:mod:`Triangulation` -- Triagulate polytope of a defined parametrization
===================================
 
.. module:: Triangulation
   :platform: Linux (Ubuntu 16.04)
   :synopsis: highlight document snippets that match a query.
.. moduleauthor:: Antonio Cao (antonio.cao@yale.edu)
 
  
Requirements::
    1.  You will need to install the numpy and scipy library to run this code.
        https://scipy.org/install.html
    2.  You will need to install SageMath to run this code:
        https://doc.sagemath.org/html/en/installation/
 
"""

import numpy as np
import scipy
import math
import os
from os.path import expanduser
from pathlib import Path
PointConfiguration.set_engine('internal')

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
    print ('pts_max and pts_min: ', pts_max, pts_min)
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
    """Compute cross product of four 4-vectors
         
        Returns a 4 vector.
         
        Args:
            v1, v2, v3, v4 are four 4-vectors
            
        Returns:
            v: the cross product of v1, v2, v3, and v4
             
        """
    #print ("input vectors: ", v1, v2, v3, v4)
    v = np.zeros((4,))
    counter = 0
    
    for i in range(4):
        mat = [v1[np.arange(len(v1))!=i].tolist(), v2[np.arange(len(v2))!=i].tolist(), v3[np.arange(len(v3))!=i].tolist()]
        mat = matrix(ZZ, mat)
        #print ('matrix: ')
        #print (mat)
        if counter == 1:
            v[i] = -1*mat.det()
            counter = 0
            #print ('neg: ', v[i])
            continue
        elif counter == 0:
            v[i] = mat.det()
            counter = 1
            #print ('pos: ', v[i])
    #print v
    mat = matrix(RR, [v1.tolist(), v2.tolist(), v3.tolist(), v4.tolist()])
    
    if mat.det() < 0:
        #print ('original: ', v)
        v = -1*v
        #print('changed: ', v)
    #print ('vector: ', v)
    return v

def Hilb(triang_list):
    """Compute the Hilbert series of a polytope.
         
        Args:
            triang_list: the triangulation of the polytope
            
        Returns:
            Series: the Hilbert series
            Triang: the triangulation (the same as input but with 1's at the end of each coordinate)
            power: the order of the t's in the Hilbert series
             
        """
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
        #print ('Hilbert: ', hilb)
        Hilb += hilb
#     print ('Hilb: ', Hilb())
    #print (Hilb(t1=t, t2=t, t3=t).series(t4, 3))
    #print ("p-q web: ", power )
    
    
    m = var('m')
    b1 = var('b1')
    b2 = var('b2')
    b3 = var('b3')
    b4 = var('b4')
    Hilb *= m^4
    
    #print ('Hilb: ', str(Hilb(t1 = (m*b1).exp(), t2 = (m*b2).exp(), t3 = (m*b3).exp(), t4 = (m*4).exp())).replace('e', 'E'))
    Series = Hilb(t1 = (m*b1).exp(), t2 = (m*b2).exp(), t3 = (m*b3).exp(), t4 = (m*4).exp()).series(m==0, 1)
    
    Series = Series.truncate()
    
    if type(power) == np.ndarray or type(power) == np.array:
        power = power.tolist()
    
    return Series, triang, power

def idx_to_pts(triang, pts):
    # Input a list of lists of indicies
    # Output a list of lists of points
    triang_new = []
    for i in range(len(triang)):
        triang_new.append([pts[j] for j in triang[i]])
    return triang_new

def lift_prism(h1, h2, h3):
    """Triangulate the polytope defined by input heights
       
       Tried Delaunay and other triangulation methods before, but none worked better than
       this simple triangulation by hand, due to constraints imposed by Physics.
         
        Args:
            h1, h2, h3: the heights of the input polytope
            
        Returns:
            prism: the point set (vertices) of the input polytope
            series: the Hilbert series
            triang: the triangulation (the same as input but with 1's at the end of each coordinate)
            power: the order of the t's in the Hilbert series
             
        """
    
    # h1 >= h2 >= h3
    
    triang = []
    power = []
    if h1 == h2 == h3:
        series = 0
        euler = 0
        prism = []
        # h1 = h2 = h3
        for h in range(h1):
            # 0 <= h < h1
            prism_1 = [[[1,0,0],[0,0,h],[0,0,h+1],[1,1,h3]]]
            prism_2 = [[[1,1,h3],[0,0,h],[0,0,h+1],[0,1,h2]]]
            prism += prism_1 + prism_2
            euler += 2
            
            series_1, triang_1, power_1 = Hilb(prism_1)
            series_2, triang_2, power_2 = Hilb(prism_2)
            series += series_1 + series_2
            triang.append(triang_1)
            triang.append(triang_2)
            power.append(power_1)
            power.append(power_2)
        return prism, series, triang, power, euler
            
    if h2 > h3:
        # h1 >= h2 > h3
        assert((h1 >= h2) and (h2 > h3))
        series = 0
        euler = 0
        prism = []
        for h in range(h3):
            # 0 <= h < h3
            prism_1 = [[[1,0,0],[0,0,h],[0,0,h+1],[1,1,h3]]]
            prism_2 = [[[0,1,h2],[0,0,h],[0,0,h+1],[1,1,h3]]]
            prism += prism_1 + prism_2
            euler += 2
            
            series_1, triang_1, power_1 = Hilb(prism_1)
            series_2, triang_2, power_2 = Hilb(prism_2)
            series += series_1 + series_2
            triang.append(triang_1)
            triang.append(triang_2)
            power.append(power_1)
            power.append(power_2)

        for h in range(h3,h2):
            # h3 <= h < h2
            prism_1 = [[[1,0,0],[0,0,h],[0,0,h+1],[1,1,h3]]]
            prism_2 = [[[0,1,h2],[0,0,h],[0,0,h+1],[1,1,h3]]]
            prism += prism_1 + prism_2
            euler += 2
            
            series_1, triang_1, power_1 = Hilb(prism_1)
            series_2, triang_2, power_2 = Hilb(prism_2)
            series += series_1 + series_2
            triang.append(triang_1)
            triang.append(triang_2)
            power.append(power_1)
            power.append(power_2)

        if h1 > h2:
            # h1 > h2 > h3
            assert(h1 > h2 and h2 > h3)
            for h in range(h2,h1):
                # h2 <= h < h1
                prism_1 = [[[1,0,0],[0,0,h],[0,0,h+1],[1,1,h3]]]
                prism_2 = [[[0,1,h2],[0,0,h],[0,0,h+1],[1,1,h3]]]
                prism += prism_1 + prism_2
                euler += 2
                
                series_1, triang_1, power_1 = Hilb(prism_1)
                series_2, triang_2, power_2 = Hilb(prism_2)
                series += series_1 + series_2
                triang.append(triang_1)
                triang.append(triang_2)
                power.append(power_1)
                power.append(power_2)
        return prism, series, triang, power, euler
        
    elif h2 < h3:
        series = 0
        euler = 0
        prism = []
        if h1 > h3:
            # h1 > h3 > h2
            assert(h1 > h3 and h3 > h2)
            for h in range(h2):
                # 0 <= h < h2
                prism_1 = [[[1,0,0],[0,0,h],[0,0,h+1],[0,1,h2]]]
                prism += prism_1
                euler += 1
                
                series_1, triang_1, power_1 = Hilb(prism_1)
                series += series_1
                triang.append(triang_1)
                power.append(power_1)
            
            # Middle prism:
            prism_1 = [[[1,0,0],[0,1,h2],[0,0,h2],[1,1,h3]]]
            prism += prism_1
            euler += 1
            
            series_1, triang_1, power_1 = Hilb(prism_1)
            series += series_1
            triang.append(triang_1)
            power.append(power_1)
            
            for h in range(h2,h3):
                # h2 <= h < h3
                prism_1 = [[[1,0,0],[0,0,h],[0,0,h+1],[1,1,h3]]]
                prism_2 = [[[0,1,h2],[0,0,h],[0,0,h+1],[1,1,h3]]]
                prism += prism_1 + prism_2
                euler += 2
                
                series_1, triang_1, power_1 = Hilb(prism_1)
                series_2, triang_2, power_2 = Hilb(prism_2)
                series += series_1 + series_2
                triang.append(triang_1)
                triang.append(triang_2)
                power.append(power_1)
                power.append(power_2)
            return prism, series, triang, power, euler
                
        else:
            # h1 = h3 > h2
            assert(h1 == h3 and h3 > h2)
            for h in range(h2):
                # 0 <= h < h2
                prism_1 = [[[1,0,0],[0,0,h],[0,0,h+1],[0,1,h2]]]
                prism += prism_1
                euler += 1
                
                series_1, triang_1, power_1 = Hilb(prism_1)
                series += series_1
                triang.append(triang_1)
                power.append(power_1)
                
            # Middle prism:
            prism_1 = [[[1,0,0],[0,1,h2],[0,0,h1],[1,1,h3]]]
            prism += prism_1
            euler += 1
            
            series_1, triang_1, power_1 = Hilb(prism_1)
            series += series_1
            triang.append(triang_1)
            power.append(power_1)
                
            for h in range(h2,h1):
                # h2 <= h < h1
                prism_1 = [[[1,0,0],[0,0,h],[0,0,h+1],[0,1,h2]]]
                prism += prism_1
                euler += 1
                
                series_1, triang_1, power_1 = Hilb(prism_1)
                series += series_1
                triang.append(triang_1)
                power.append(power_1)
            return prism, series, triang, power, euler
                
    else:
        assert(h2 == h3)
        assert(h1 > h2)
        # h1 > h2 = h3
        prism = []
        series = 0
        euler = 0
        for h in range(h2):
            # 0 <= h < h2
            prism_1 = [[[1,0,0],[0,0,h],[0,0,h+1],[1,1,h3]]]
            prism_2 = [[[0,1,h2],[0,0,h],[0,0,h+1],[1,1,h3]]]
            prism += prism_1 + prism_2
            euler += 2

            series_1, triang_1, power_1 = Hilb(prism_1)
            series_2, triang_2, power_2 = Hilb(prism_2)
            series += series_1 + series_2
            triang.append(triang_1)
            triang.append(triang_2)
            power.append(power_1)
            power.append(power_2)
        
        for h in range(h2,h1):
            # h2 <= h < h1
            prism_1 = [[[1,0,0],[0,0,h],[0,0,h+1],[1,1,h3]]]
            prism_2 = [[[0,1,h2],[0,0,h],[0,0,h+1],[1,1,h3]]]
            prism += prism_1 + prism_2
            euler += 2

            series_1, triang_1, power_1 = Hilb(prism_1)
            series_2, triang_2, power_2 = Hilb(prism_2)
            series += series_1 + series_2
            triang.append(triang_1)
            triang.append(triang_2)
            power.append(power_1)
            power.append(power_2)
        
    return prism, series, triang, power, euler

def prism_plot(prism):
    """Plot the input prism
         
        Args:
            prism: the point set (i.e. vertices) of the input polytope
            
        Returns:
            plot: the SageMath plot of the input polytope
            (Note: this require an interactive session)
             
        """
    count = 0
    plot = 0
    for p in prism:
        poly = Polyhedron(vertices=p)
        if count == 0:
            count += 1
            plot = poly.plot()
        else:
            plot += poly.plot()
    
#     print(poly_list.vertices())
    return plot

def write_file(fname, a):
    p = Path(fname)
    test_file = p.open('ab')
    np.save(test_file, a)
    test_file.close()

def load_file(fname):
    p = Path(fname)
    out_list = []
    with p.open('rb') as f:
        fsz = os.fstat(f.fileno()).st_size
        out = np.load(f)
        out_list.append(out)
        while f.tell() < fsz:
            out_list.append(np.load(f))
    return out_list
   
def generate_series(h_max):
    """Generate triangulation and hilbert series for all polytope of max hight less than h_max
        
        The return file is saved at '~/data/' with name 'series_0_%d.npy'%(h_max-1)
         
        Args:
            h_max: the (open) upper bound of the maximum height of all the generated polytopes
            
        Note: each entry of the output data file contains
             [h1,h2,h3]: the heights
             prism: polytope point set
             series: Hilbert series
             triang: triangulation of the polytope
             power: order of the t's in the hilbert series
        """
    home = expanduser("~/Calabi_Yau")
    file_name='series'+'_0_'+str(int(h_max-1))+'.npy'
    path_name = home+'/data/'+file_name
    for h1 in range(1,h_max):
        for h2 in range(h1+1):
            for h3 in range(h2+1):
                print(h1,h2,h3)
                prism, series, triang, power, euler = lift_prism(h1,h2,h3)
                data = [str([h1,h2,h3]), str(prism), str(series), str(triang), str(power), str(euler)]
                write_file(path_name, data)

def load_series(file_name):
    home = expanduser("~/Calabi_Yau")
    path_name = home+'/data/'+file_name
    return load_file(path_name)



generate_series(50)
# data_list = load_series('series_0_2.npy')

# for data in data_list:
#     heights = eval(data[0])
#     series = data[2]
#     euler = eval(data[-1])
#     print heights, euler, series
    
