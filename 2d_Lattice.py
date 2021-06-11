import numpy as np
import itertools as it
import copy
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import pandas as pd
import csv

'error tolerance for reduction and lattice vectors'
err = 10**(-14)

'-------------------------------------------------------------'
'UTILITY FUNCTIONS
'-------------------------------------------------------------'
'Swap two elements of a list given their indices'
def swapper(list, i, j):

    if i > len(list) or j > len (list):
        print('Index error!')

    else:
        list[j], list[i] = list[i], list[j]

    return list

'Given a list, an index and a value, rebalance will'
'substract the value from the list element at the input index'
'and add the value to all other elements in the list'
def rebalance(lst, indx, val):

    for i in range (0, len(lst)):
        if i == indx:
            lst[i] -= val
        else:
            lst[i] += val

    return lst

'Standard Distance Calculation'
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

'Return Cartesian Vectors given length/angle parameters'
def angcalc_2d(a=1, b=1, alpha):
    
    return [[a,0], [b*np.cos(alpha), b*np.sin(alpha)]
    


'-------------------------------------------------------------'
'2D Lattice Object'
'-------------------------------------------------------------'
class Lattice_2d:

    def __init__(self,lattice):
        self.lattice = lattice
        self.x = np.array(self.lattice[0])
        self.y = np.array(self.lattice[1])
        self.inner = np.dot(self.x, self.y)
        self.angle = np.arccos(self.inner)
        self.mat = np.transpose(np.array(self.lattice))
        self.vol = np.linalg.det(self.mat)

    'Return a (not necessarily obtuse) superbase')
    def msb(self):

        return[-np.add(self.x, self.y), self.x, self.y]

'-------------------------------------------------------------'
'SUPERBASE OBJECT'
'-------------------------------------------------------------'
class Obsuper_2d:

    def __init__(self,basis):
        self.basis = basis
        self.vo = np.array(basis[0])
        self.vl = np.array(basis[1])
        self.vt = np.array(basis[2])
        self.p_ol = np.dot(self.vo, self.vl)
        self.p_ot = np.dot(self.vo, self.vt)
        self.p_lt = np.dot(self.vl, self.vt)
        self.sellings_t = [-self.p_ol, -self.p_ot, -self.p_lt]

    def red2d(self):

        fin_coform = copy.copy(self.conorms)
        redcount = 0
        corrcount = 0

        while min(fin_coform) < 0:
            for i in range(0,3):
                if fin_coform[i] < -err:
                    rebalance(fin_coform, i, 2*fin_coform[i])
                    redcount += 1

                elif -err < fin_coform[i] < err:
                    fin_coform[i] = 0
                    corrcount += 1

        print('number of corrections = '+ repr(corrcount))
        print('number of reduction steps = '+ repr(redcount))

        print(fin_coform)
        return fin_coform
            
'-------------------------------------------------------------'
'ISOMETRIC DISTANCE CALCULATION'
'-------------------------------------------------------------'
def isodist_2d(lat_1, lat_2, dtype = 2, reflect = False):

    coform_1 = Obsuper_2d(lat_1.msb()).red2d()
    coform_2 = Obsuper_2d(lat_2.msb()).red2d()

    if reflect:
        perlist_1 = [[coform_1[j-i] for i in range (len(coform_1))] for j in range (len(coform_1))]
        perlist_2 = [[coform_2[j-i] for i in range (len(coform_2))] for j in range (len(coform_2))]

    else:
        perlist_1 = list(it.permutations(coform_1))
        perlist_2 = list(it.permutations(coform_2))

    pairdists = []
    for i in perlist_1:
        for j in perlist_2:
            pairdists.append(dst(i,j, dtype))

    return min(pairdists)
