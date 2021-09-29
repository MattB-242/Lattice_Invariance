import numpy as np
import itertools as it
import copy
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import pandas as pd
import csv
from pymatgen import core


'-----------------------------------------------------------------'
'GENERAL UTILITY FUNCTIONS'
'-----------------------------------------------------------------'
'Check if a value is within tolerance of another value and replace if necessary'
def tolcheck (a, val, tol = 10**-14):
    if np.abs(a-val) < tol:
        return val
    else:
        return a


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

'Calculate a reduced lattice object from parameters'
def lat_param(a,b,c,alpha,beta,gamma):

    orlat = core.lattice.Lattice.from_parameters(a,b,c,alpha,beta,gamma)
    niggli_lattice = orlat.get_niggli_reduced_lattice()
    return [i for i in list(niggli_lattice.matrix)]

'perform a single reduction step on an arbitrary 2x3 array'
def one_red_step(array):

    j = array[0]
    k = array[1]
    l = j+k
    if min(j+k) < 0:
        epsilon = -min(j+k)
        if l.index(-epsilon) == 0:
            return [[epsilon, l[1]-epsilon, l[4]-epsilon], [l[3]+epsilon, l[2]-epsilon, l[5]-epsilon]]
        elif l.index(-epsilon) == 1:
            return [[l[0]-epsilon, epsilon, l[3]-epsilon], [l[2]-epsilon, l[4]+epsilon, l[5]-epsilon]]
        elif l.index(-epsilon) == 2:
            return [[l[0]-epsilon, l[3]-epsilon, epsilon], [l[1]-epsilon, l[4]-epsilon, l[5]+epsilon]]
        elif l.index(-epsilon) == 3:
            return [[l[0]+epsilon, l[1]-epsilon, l[2]-epsilon], [epsilon, l[5]-epsilon, l[4]+epsilon]]
        elif l.index(-epsilon) == 4:
            return [[l[0]-epsilon, l[1]+epsilon, l[2]-epsilon], [l[5]-epsilon, epsilon, l[3]+epsilon]]
        else:
            return [[l[0]-epsilon, l[1]-epsilon, l[2]+epsilon], [l[4]-epsilon, l[3]-epsilon, epsilon]]

    else:
        return [l,k]

'Generally rearrange a 2x3 array into coform-canonical order'
def canon_coform(array, err = 10**-14):

    #place smallest value at [0][0]
    cof_red = list(zip(array[0],array[1]))
    cof_sort = sorted(cof_red, key = lambda val : val[0])
    sk = [i[0] for i in cof_sort]
    sl = [i[1] for i in cof_sort]
    right_col = [sk[1], sk[2], sl[1], sl[2]]

    #place next smallest value at [0][1]
    maxind = right_col.index(min(right_col))
    if maxind == 0:
        return [[tolcheck(i,0,err) for i in sk],[tolcheck(i,0,err) for i in sl]]
    elif maxind == 1:
        sk_mod = [sk[0], sk[2], sk[1]]
        sl_mod = [sl[0], sl[2], sl[1]]
    elif maxind == 2:
        sk_mod = [sk[0], sl[1], sl[2]]
        sl_mod= [sl[0], sk[1], sk[2]]
        return [[tolcheck(i,0,err) for i in sk_mod],[tolcheck(i,0,err) for i in sl_mod]]
    else:
        sk_mod = [sk[0], sl[2], sl[1]]
        sl_mod= [sl[0], sk[2], sk[1]]
        return [[tolcheck(i,0,err) for i in sk_mod],[tolcheck(i,0,err) for i in sl_mod]]

'Calculate barycentric co-ordinates in the canonical triangle from a list of three positive reals'
def baryproj(l):

    s = sum(l)
    return [(l[2]-l[1])/(2*s), l[0]/s]

'Calculate orthogonal projection along triangular symmetric unit vectors from a list of three positive reals'
def orthproj(l):

    s = sum(l)/3

    return [(l[1]-l[0])*np.sqrt(3)/2, (l[2]-s)-(l[1]-s)/2-(l[0]-s)/2]

'-----------------------------------------------------------------'
'COFORM, VOFORM and LATTICE OBJECTS'
'-----------------------------------------------------------------'


'3D LATTICE OBJECT'
'--------------------------------------------------------------'

'Takes three basis vectors as input'

class Lattice_3d:

    def __init__(self,lat):
        self.lat = lat
        self.x = np.array(lat[0])
        self.y = np.array(lat[1])
        self.z = np.array(lat[2])
        self.allvects = [self.x, self.y, self.z]
        self.a = np.linalg.norm(lat[0])
        self.b = np.linalg.norm(lat[1])
        self.c = np.linalg.norm(lat[2])
        self.al = np.arccos(np.dot(self.y, self.z)/(self.b*self.c))
        self.bet = np.arccos(np.dot(self.x, self.z)/(self.a*self.c))
        self.gam = np.arccos(np.dot(self.x, self.y)/(self.b*self.c))
        self.mat = np.transpose(np.array(self.lat))
        self.gram = np.matmul(np.array(self.lat), self.mat)

    'Niggli reduce the input lattice basis'
    def latnig(self):

        return Lattice_3d(lat_param(self.a, self.b, self.c, self.al, self.bet, self.gam))

    'generates a putative coform out of the initial cartesian vectors'
    def makecof(self):

        xy = np.add(self.x, self.y)

        l = [-np.add(xy, self.z), self.x, self.y, self.z]

        return Coform([[-np.dot(l[2], l[3]), -np.dot(l[1], l[3]), -np.dot(l[1], l[2])],
                [-np.dot(l[0], l[1]), -np.dot(l[0], l[2]), -np.dot(l[0], l[3])]])


    'Get the lattice coform in standard order directly'
    def get_lat_coform(self):

         return Coform(self.latnig().makecof()).red_cof()

    def get_lat_rootform(self):

        return Coform(self.get_lat_coform().cf_root())

'----------------------------------------------------------------'
'COFORM OBJECT: input format is a 2 x 3 array'
'----------------------------------------------------------------'
class Coform:

    def __init__ (self, conorms):
        self.conorms = list(conorms)
        self.toprow = self.conorms[0]
        self.bottomrow = self.conorms[1]
        self.asvect = self.toprow + self.bottomrow
        self.p_23 = self.conorms[0][0]
        self.p_13 = self.conorms[0][1]
        self.p_12 = self.conorms[0][2]
        self.p_01 = self.conorms[1][0]
        self.p_02 = self.conorms[1][1]
        self.p_03 = self.conorms[1][2]

    def get_voform(self):
        v_0 = sum(self.bottomrow)
        v_1 = self.p_01 + self.p_12 + self.p_13
        v_2 = self.p_02 + self.p_12 + self.p_23
        v_3 = self.p_03 + self.p_13 + self.p_23
        v_12 = self.p_01 + self.p_13 + self.p_02 + self.p_23
        v_13 = self.p_01+ self.p_12 + self.p_03 + self.p_23
        v_23 = self.p_02 + self.p_12 + self.p_03 + self.p_13

        return Voform([[v_0], [v_1, v_2, v_3], [v_23, v_13, v_12]])

    def red_cof(self):
        putative = self.conorms

        while min(putative[0] + putative[1]) < 0:
            putative = one_red_step(putative)

        return Coform(canon_coform(putative))


    'Puts coform in canonical form'
    def cf_canonical(self):

            return Coform(canon_coform(self.conorms))

    'Returns root_form'
    def cf_root(self):

        return [[np.sqrt(i) for i in self.toprow], [np.sqrt(i) for i in self.bottomrow]]

    'General barycentric coform plot co-ordinates'
    def cf_root_bary(self):

        return baryc(self.cf_root()[0]) + baryc(self.cf_root[1])

    'Outputs list of all coform permutations'
    def cf_permute(self):

            s = self.conorms

            s_dict = {(2,3):self.p_23, (1,3):self.p_13, (1,2):self.p_12,
                      (0,1):self.p_01, (0,2):self.p_02, (0,3):self.p_03}

            base_list = [0,1,2,3]
            perms_list = [list(i) for i in it.permutations(base_list)]


            all_perms = []
            for j in perms_list:
                    perm_dict = {}
                    for k in s_dict.keys():
                        #permute indices according to j and find the right key
                        ind_list = [base_list[j.index(k[0])], base_list[j.index(k[1])]]
                        ind_tup = tuple(sorted(ind_list))
                        perm_dict[k] = s_dict[ind_tup]
                    all_perms.append([perm_dict[(2,3)], perm_dict[(1,3)], perm_dict[(1,2)],
                                      perm_dict[(0,1)], perm_dict[(0,2)], perm_dict[(0,3)]])


            return all_perms

'----------------------------------------------------------------'
'VOFORM OBJECT = Input format is a list containing one element , followed by a 2x3 array'
'----------------------------------------------------------------'
class Voform:

    def __init__(self, vonorms):
        self.vonorms = list(vonorms)
        self.volist = self.vonorms[0] + self.vonorms[1] + self.vonorms[2]
        self.lat_lengths = [np.sqrt(i) for i in self.vonorms[1]]
        self.sum_lengths = [np.sqrt(i) for i in self.vonorms[2]]
        self.v_0 = self.vonorms[0][0]
        self.v_1 = self.vonorms[1][0]
        self.v_2 = self.vonorms[1][1]
        self.v_3 = self.vonorms[1][2]
        self.v_23 = self.vonorms[2][0]
        self.v_13 = self.vonorms[2][1]
        self.v_12 = self.vonorms[2][2]

    'Calculates coform from voform'
    def get_coform(self):
        p_01 = (self.v_0 + self.v_1 - self.v_23)/2
        p_02 = (self.v_0 + self.v_2 - self.v_13)/2
        p_03 = (self.v_0 + self.v_3 - self.v_12)/2
        p_12 = (self.v_1 + self.v_2 - self.v_12)/2
        p_13 = (self.v_1 + self.v_3 - self.v_13)/2
        p_23 = (self.v_2 + self.v_3 - self.v_23)/2

'-----------------------------------------------------------'
'COFORM-VOFORM PAIR'
'Entered as a 2x3 array followed by a list of six numbers'
'-----------------------------------------------------------'
class CofVof:

    def __init__(self, cofvof):
        self.cofvof = list(cofvof)
        self.coform = self.cofvof[0]
        self.voform = self.cofvof[1]


    def check(self):
        if self.coform.get_voform() != self.voform.vonorms:
            print("Coform and voform are not from the same lattice!")
        else:
            print("Coform-Voform pair checks out!")

    def get_angles(self):

        return[np.rad2deg(np.arccos(-self.coform.p_23/(self.get_lat_lens()[1]*self.get_lat_lens()[2]))),
               np.rad2deg(np.arccos(-self.coform.p_13/(self.get_lat_lens()[0]*self.get_lat_lens()[2]))),
               np.rad2deg(np.arccos(-self.coform.p_12/(self.get_lat_lens()[0]*self.get_lat_lens()[1])))]

    def get_lattice_vectors(self):
        parameters = self.voform.lat_lengths + self.get_angles()
        #print('input parameters are ' + repr(parameters))
        unreduced_lattice = core.lattice.Lattice.from_parameters(*parameters)
        return [list(i) for i in unreduced_lattice.matrix]

    def get_lat_lens(self):

        return self.voform.lat_lengths

    def red_cv_one(self):

        if min(self.coform.asvect) < 0:
            ep = min(self.coform.asvect)
            ind = self.coform.asvect.index(ep)

            if ind == 1:
                newcof = Coform(one_red_step(self.coform.conorms))
                newvof = Voform([[self.voform.v_23], [self.voform.v_1, self.voform.v_12, self.voform.v_3],
                        [self.voform.v_0, self.voform.v_13 +(4*ep), self.voform.v_2]])

            elif ind == 0:
                newcof = Coform(one_red_step(self.coform.conorms))
                newvof = Voform([[self.voform.v_13], [self.voform.v_12, self.voform.v_2, self.voform.v_3],
                        [self.voform.v_23+(4*ep), self.voform.v_0, self.voform.v_1]])

            elif ind == 2:
                newcof = Coform(one_red_step(self.coform.conorms))
                newvof = Voform([[self.voform.v_23], [self.voform.v_1, self.voform.v_2, self.voform.v_13],
                        [self.voform.v_0, self.voform.v_3, self.voform.v_12+ (4*ep)]])

            elif ind == 3:
                newcof = Coform(one_red_step(self.coform.conorms))
                newvof = Voform([[self.voform.v_0], [self.voform.v_1, self.voform.v_13, self.voform.v_12],
                        [self.voform.v_23+(4*ep), self.voform.v_2, self.voform.v_3]])


            elif ind == 4:
                newcof = Coform(one_red_step(self.coform.conorms))
                newvof = Voform([[self.voform.v_0], [self.voform.v_23, self.voform.v_2, self.voform.v_12],
                        [self.voform.v_1, self.voform.v_2 + (4*ep), self.voform.v_3]])

            elif ind == 5:
                newcof = Coform(one_red_step(self.coform.conorms))
                newvof = Voform([[self.voform.v_0], [self.voform.v_23, self.voform.v_13, self.voform.v_3],
                        [self.voform.v_1, self.voform.v_2, self.voform.v_12+(4*ep)]])

        else:
                newcof = self.coform
                newvof = self.voform


        return CofVof([newcof, newvof])

    def full_reduction(self):

        while min(self.coform.asvect) < 0:

            self = self.red_cv_one()

        return self


'--------------------------------------------------------------'
'LATTICE DISTANCE CALCULATION FUNCTIONS'
'--------------------------------------------------------------'
'Takes two lattices as input. Use dtype = 1/2 for L_1/L_2 distance, dtype = 0 for L_inf'

def lat_cfdist(lat_1, lat_2, dtype = 2):

    oscof_1 = lat_1.get_lat_rootform()
    oscof_2 = lat_2.get_lat_rootform()


    osperms_2 = oscof_2.cf_permute()

    alldists = []
    for i in osperms_2:
        alldists.append(dst(oscof_1.toprow + oscof_1.bottomrow,i,type = dtype))

    return min(alldists)

'calculate distance between two input coforms'
def cfdist(cof_1, cof_2, dtype = 2):

    dist_1 = cof_1.toprow + cof_1.bottomrow
    dist_2 = cof_2.cf_permute()

    for i in dist_2:
        distlist.append(dst(dist_1, i, type = dtype))

    return(min(distlist))
