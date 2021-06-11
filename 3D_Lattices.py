'-----------------------------------------------------------------'
'IMPORTS'
'-----------------------------------------------------------------'

import numpy as np
import itertools as it
import copy
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import pandas as pd
import csv

'Error tolerance for lattice vectors and coform reduction'
err = 10**(-14)

'-----------------------------------------------------------------'
'GENERAL UTILITY FUNCTIONS
'-----------------------------------------------------------------'

'Swap two elements of a list given their indices'
def swapper(list, i, j):

    if i > len(list) or j > len (list):
        print('Index error!')

    else:
        list[j], list[i] = list[i], list[j]

    return list

'Given a list, an index and a value,'
'subtract the value from the list element at the input index'
'and add the value to all other elements in the list'
def rebalance(lst, indx, val):

    for i in range (0, len(lst)):
        if i == indx:
            lst[i] -= val
        else:
            lst[i] += val

    return lst


'General distance function for vectors. Takes type input:'
'1 = L1, 2= L2, 0 = L_inf'
def dst(a,b,type):

    if len(a) != len(b):
        print('Both lists must be the same length')

        return 0

    elif type == 1:
        tot = 0
        for i in range(0, len(a)):
            tot += np.abs(a[i] - b[i])

        return tot

    elif type == 2:
        return np.linalg.norm(np.array(b)-np.array(a))

    else:
        dlist = []
        for i in range (0, len(a)):
            dlist.append(np.abs(a[i] - b[i]))

        return max(dlist)
      
'----------------------------------------------------------------------------'
'HARDCODED BRAVAIS LATTICE TYPES
'----------------------------------------------------------------------------'

'Generates Lattice with input a, b, c, alpha, beta, gamma parameters'
'Allowed Lattices are Cubic, Hexagonal, Rhombohedral, Tetragonal, Orthorhombic and Triclinic'
'Allowed Lattice Types are P = Primitive, C = Base Centered, F = Face Centered, I = Body Centered'
'If a disallowed type is used type defaluts to P. If a lattice is not recognised defaults to cubic'
'All parameters must be provided. Cubic uses only the value of a, tetragonal only a and c'
'orthorhombic uses the beta angle value and defaults to 90 otherwise. 

def cartcalc_3d(lattice, a, b, c, alpha, beta, gamma, type = 'P'):
    if lattice == 'Cubic':

        if type == 'P':
            return [[a,0,0],[0,a,0], [0,0,a]]
        elif type == 'F':
            return [[0,a/2,a/2], [a/2, 0, a/2], [a/2, a/2,0]]
        elif type == 'I':
            return [[-a/2,a/2,a/2], [a/2, -a/2, a/2], [a/2, a/2,-a/2]]
        else:
            print('Lattice Subtype Not Recognised, Defaulting to Cubic')
            return [[a, 0, 0], [0,a,0], [0,0,a]]

    elif lattice == 'Hexagonal':
        return [[a/2, - (a*np.sqrt(3))/2, 0], [a/2, (a*np.sqrt(3))/2, 0], [0,0,c]]

    elif lattice == 'Rhombohedral':
        return [[a/2, -a/(2*np.sqrt(3)), c/3], [0, a/np.sqrt(3), c/3], [-a/2, -a/(2*np.sqrt(3)), c/3]]

    elif lattice == 'Tetragonal':

        if a == c:
            print('This is a cubic lattice!')

        if type == 'P':
            return [[a,0,0], [0,a,0], [c,0,0]]
        elif type == 'I':
            return [[-a/2, a/2, c/2], [a/2, -a/2, c/2], [a/2, a/2, -c/2]]
        else:
            print('Type not recognised, defaulting to primitive')

    elif lattice == 'Orthorhombic':

        if type == 'P':
            return [[a,0,0], [0,b,0], [0,0,c]]
        elif type == 'C':
            return [[a/2, -b/2, 0], [a/2, b/2, 0], [0,0,c]]
        elif type == 'I':
            return [[-a/2, b/2, c/2], [a/2, -b/2, c/2], [a/2, b/2, c/2]]
        elif type == 'F':
            return[[0,b/2,c/2],[a/2,0,c/2], [a/2,b/2,0]]
        else:
            print('Type not recognised, defaulting to primitive')
            return [[a,0,0], [0,b,0], [0,0,c]]

    elif lattice == 'Monoclinic':

        if beta == 90:
            print('Lattice is Orthorhombic!')

        if type == 'P':
            return[[a,0,0], [0,b,0], [c*np.cos(np.deg2rad(beta)), 0, c*np.sin(np.deg2rad(beta))]]
        elif type == 'C':
            return [[a/2,-b/2, 0], [a/2, b/2, 0], [c*np.cos(np.deg2rad(beta)), 0, c*np.sin(np.deg2rad(beta))]]
        else:
            print('Type Not Recognised, Reverting to Primitive')
            return[[a,0,0], [0,b,0], [c*np.cos(np.deg2rad(beta)), 0, c*np.sin(np.deg2rad(beta))]]

    elif lattice == 'Triclinic':

        c_x = c*np.cos(np.deg2rad(beta))
        c_y = c*((np.cos(np.deg2rad(alpha))-(np.cos(deg2rad(beta)*np.cos(deg2rad(gamma)))))/np.sin(np.deg2rad(gamma)))
        c_z = np.sqrt(c**2 -c_x**2 = c_y**2)

        return [[a,0,0], [b*np.cos(np.deg2rad(gamma)), b*np.sin(np.deg2rad(gamma)), 0], [c_x, c_y, c_z]]

    else:
        print('Lattice Category Not Recognised, Defaulting to Primitive Cubic')
        return[[a,0,0], [0,a,0], [0,0,a]]
 
'--------------------------------------------------------------'
'3D LATTICE OBJECT'
'--------------------------------------------------------------'

'Cartesian Vectors Required - use cartcalc_3d if you have length and angle parameters'

class Lattice_3d:

    def __init__(self,lattice):
        self.lattice = lattice
        self.x = np.array(lattice[0])
        self.y = np.array(lattice[1])
        self.z = np.array(lattice[2])
        self.mat = np.transpose(np.array(self.lattice))
        self.gram = np.matmul(np.array(self.lattice), self.mat)

    'generates a superbase out of the initial cartesian vectors
    def msb3d(self):

        xy = np.add(self.x, self.y)

        return[-np.add(xy, self.z), self.x, self.y, self.z]
      
 '--------------------------------------------------------------'
'SUPERBASE OBJECT'
'--------------------------------------------------------------'

'Takes output of msb3d lattice object function'

class Obsuper_3d:

    def __init__(self,basis):
        self.basis = basis
        self.vo = np.array(self.basis[0])
        self.vl = np.array(self.basis[1])
        self.vt = np.array(self.basis[2])
        self.vh = np.array(self.basis[3])
        self.vol = np.add(self.vo, self.vl)
        self.vot = np.add(self.vo, self.vt)
        self.voh = np.add(self.vo, self.vh)
        self.p_ol = np.dot(self.vo, self.vl)
        self.p_ot = np.dot(self.vo, self.vt)
        self.p_oh = np.dot(self.vo, self.vh)
        self.p_lt = np.dot(self.vl, self.vt)
        self.p_lh = np.dot(self.vl, self.vh)
        self.p_th = np.dot(self.vt, self.vh)
        self.put_vonorms = [np.dot(self.vo,self.vo), np.dot(self.vl,self.vl), np.dot(self.vt,self.vt),
                        np.dot(self.vh, self.vh), np.dot(self.vol,self.vol), np.dot(self.vot,self.vot),
                        np.dot(self.voh,self.voh)]
        self.put_conorms = [-self.p_ol, -self.p_ot, - self.p_oh, -self.p_lt, -self.p_lh, -self.p_th]



    'Reduction method as outlined in Conway and Sloane'
    'Returns two triples as input for distance measurement'
    def reduce_3d(self):

        cent = 0
        l = [self.p_lt, self.p_lh, self.p_th]
        k = [self.p_ol, self.p_ot, self.p_oh]


        while max(l) > 0 or cent > 0:

            if cent > 0:
                for i in range(0,3):
                    if l[i] == 0:
                        rebalance(l, i, cent)
                        rebalance (k, 2-i, cent)
                        cent -= cent

            elif cent < 0:
                j = max(l)
                cent += j
                l = [i - j for i in l]
                k = [i + j for i in k]

            else:
                j = max(l)
                ind = l.index(j)
                l = rebalance(l, ind, j)
                k = rebalance(k, 2-ind, j)
                cent -= j

        return [[-k[0], -k[1], -k[2]], [-l[2], -l[1], -l[0]]]

    'Outputs a list of all coform permutations'
    def allperms_3d(self, iso = False):

            permlist = []
            swapind = [[0,1], [0,2], [1,2]]

            a = self.reduce_3d()[0]
            b = self.reduce_3d()[1]

            if iso:
                pa = [[a[j-i] for i in range (len(a))] for j in range (len(a))]
                pb = [[b[j-i] for i in range (len(b))] for j in range (len(b))]

            else:
                pa = [list(i) for i in list(it.permutations(a))]
                pb = [list(j) for j in list(it.permutations(b))]

            permlist = [pa[i] + pb[i] for i in range (0,len(pa))]


            swaplist = []
            for i in permlist:
                swp = copy.copy(i)
                for j in swapind:
                    for k in j:
                        swapper(swp, k, k+3)
                    swaplist.append(swp)

            return permlist + swaplist
          

 '--------------------------------------------------------------'
'ISOMETRIC DISTANCE CALCULATION FUNCTION'
'--------------------------------------------------------------'
'Takes two lattices as input. Use dtype = 1/2 for L_1/L_2 distance, dtype = 0 for L_inf'
'Set iso = True to suppress coform permutations induced by lattice reflection'
def lattice_3d_distance(lat_1, lat_2, dtype = 2, iso = False):

    osbase_1 = Obsuper_3d(lat_1.msb3d())
    osbase_2 = Obsuper_3d(lat_2.msb3d())

    osperms_1 = osbase_1.allperms_3d(iso = iso)
    osperms_2 = osbase_2.allperms_3d(iso = iso)

    alldists = []
    for i in osperms_1:
        for j in osperms_2:
            alldists.append(dst(i,j, type = dtype))

    return min(alldists)
