import numpy as np
import pandas as pd
import csv
import copy
from math import atan2
from math import degrees

infile = '/Users/mattbright/Desktop/2D_Data/all_CSD_nodisord_cent_870630_clean_coforms.csv'
outpath = '/Users/mattbright/Desktop/2D_Data'

tol = 10**-8

'-----------------------------------------------------'
'UTILITY FUNCTIONS'
'-----------------------------------------------------'
'General L_q distance calculator for vectors in R_2'
def distgen(q, v_1, v_2):

    if q == 0:
        return max(np.abs(v_1[0] - v_2[0]), np.abs(v_1[1] - v_2[1]))
    else:
        return ((np.abs(v_1[0] - v_2[0])**q) + (np.abs(v_1[1] - v_2[1])**q))**(1/q)

'A single reduction step on a three member list'
def redstep(l, index):
    newlist = []
    eps = l[index]
    for i in l:
        if l.index(i) == index:
            newlist.append(-eps)
        else:
            newlist.append(i+(2*eps))

    return newlist

'Cycle any list so that its least member is first'
def mincyc(l):

    top = l.index(min(l))

    return l[top:] + l[:top]

def roundlist(l, r=2):

    rl = []
    for i in l:
        rl.append(round(i,r))

    return rl

'Calculate M_inf quantity (see paper, Lemma C.2)'
def minf(a,b,c,d):
    return max([np.abs(a-b), np.abs(c-d), (np.abs(a+b-c+d))/2])

'Make a 2D Lattice that from length and angle parameters (in degrees)'
def makelat (a, b, t):
    trad = np.deg2rad(t)

    return Lat2d([a, 0], [b*np.cos(trad), b*np.sin(trad)])

def countsame(l):
    count = [l.count(i) for i in l]
    return max(count)

def sb_sign(veclist):

    u_1 = veclist[1]/np.linalg.norm(veclist[1])
    u_0 = veclist[0]/np.linalg.norm(veclist[0])
    u_2 = veclist[2]/np.linalg.norm(veclist[2])

    ang_1 = degrees(atan2(u_1[1], u_1[0])) % 360
    ang_0 = degrees(atan2(u_0[1], u_0[0])) % 360
    ang_2 = degrees(atan2(u_2[1], u_2[0])) % 360

    ang_12 = (ang_2 - ang_1) % 360
    ang_10 = (ang_0 - ang_1) % 360

    if ang_12 < ang_10:
        sgn = 1
    else:
        sgn = -1

    return ([v_0, v_1, v_2], sgn)

'----------------------------------------------------'
'LATTICE CLASS IN 2D WITH ROOT FORM GENERATOR'
'----------------------------------------------------'
class Lat2d:

    def __init__(self, vec_1, vec_2):
        self.x = np.array(vec_1)
        self.y = np.array(vec_2)
        self.ob = -(self.x + self.y)
        self.xlen = np.linalg.norm(vec_1)
        self.ylen = np.linalg.norm(vec_2)
        self.oblen = np.linalg.norm(self.ob)
        self.inner = np.dot(vec_1, vec_2)
        self.angle = np.rad2deg(np.arccos(self.inner/(self.xlen * self.ylen)))

    def make_obsb(self):

        v_1 = self.x
        v_2 = self.y
        v_0 = self.ob


        inners = [-np.dot(v_1, v_2), -np.dot(v_0, v_1), -np.dot(v_0, v_2)]
        stepcount = 0

        while min(inners) < 0:
            if inners[0] < 0:
                v_0 = v_1 - v_2
                v_1 = -v_1
            elif inners[1] < 0:
                v_2 = v_0 - v_1
                v_0 = -v_0
            else:
                v_1 = v_0 - v_2
                v_0 = -v_0
            inners = [-np.dot(v_1, v_2), -np.dot(v_0, v_1), -np.dot(v_0, v_2)]

        return [v_0, v_1, v_2]

    def sb_sign(self):

        veclist = self.make_obsb()

        u_1 = veclist[1]/np.linalg.norm(veclist[1])
        u_0 = veclist[0]/np.linalg.norm(veclist[0])
        u_2 = veclist[2]/np.linalg.norm(veclist[2])

        ang_1 = degrees(atan2(u_1[1], u_1[0])) % 360
        ang_0 = degrees(atan2(u_0[1], u_0[0])) % 360
        ang_2 = degrees(atan2(u_2[1], u_2[0])) % 360

        ang_12 = (ang_2 - ang_1) % 360
        ang_10 = (ang_0 - ang_1) % 360

        if ang_12 < ang_10:
            sgn = 1
        else:
            sgn = -1

        return sgn

    def make_cf(self):

        obsb = self.make_obsb()
        return [-np.dot(obsb[1], obsb[2]), -np.dot(obsb[0], obsb[1]), -np.dot(obsb[0], obsb[2])]


    def lattice_sign(self):

        cf = mincyc(self.make_cf())
        sgn = self.sb_sign()

        if cf[1] > cf[2]:
            sgn = -sgn
        else:
            sgn = sgn

        return sgn

    def make_rf(self):

        cf = self.make_cf()
        sgn = self.lattice_sign()

        rf = sorted([np.sqrt(i) for i in cf])

        if np.abs(rf[0])<tol or (np.abs(rf[0] - rf[1]) < tol or np.abs(rf[1] - rf[2])<tol):
                return RF2_signed(rf, 0)
        else:
                return RF2_signed(rf, sgn)



'----------------------------------------------------'
'2D ORIENTED ROOT FORM CLASS'
'----------------------------------------------------'

class RF2_signed:

    def __init__(self, vec, sign):
        self.vec = vec
        self.r_12 = vec[0]
        self.r_01 = vec[1]
        self.r_02 = vec[2]
        self.sign = sign

        'Throws a warning if the root form is not ordered'
        if sorted(self.vec) != self.vec:
            print('Warning! Root form is not ordered.')

    'Create correctly ordered unoriented root form'
    def rightsign(self):
        if self.sign == -1:
            return [self.vec[0], self.vec[2], self.vec[1]]
        else:
            return self.vec

    'Calculate Chirality (L_inf and L_2)'
    def rf_grpchir(self, pgroup = 2, dtype = 0):
        if pgroup == 2:
            if dtype == 0:
                return min([self.r_12, (self.r_01 - self.r_12)/2, (self.r_02 - self.r_01)/2])
            elif dtype == 2:
                return min([self.r_12, (self.r_01 - self.r_12)/np.sqrt(2), (self.r_02 - self.r_01)/np.sqrt(2)])
            else:
                print('I can only calculate L_2 and L_inf distances')
                return 0
        elif pgroup == 4:
            if dtype == 0:
                return max([2*self.r_12, (self.r_01-self.r_02)/2])
            elif dtype == 2:
                return (np.sqrt((4*self.r_12) +(self.r_01 + self.r_02)**2))/2
                print('I can only calculate L_2 and L_inf distances')
                return 0
        elif pgroup == 6:
            if dtype == 0:
                return self.r_02 - self.r_12
            elif dtype == 2:
                return np.sqrt(2/3 *(self.r_12**2 + self.r_01**2 + self.r_02**2 -(self.r_01*self.r_12) - (self.r_12*self.r_02) - (self.r_01*self.r_02)))
                print('I can only calculate L_2 and L_inf distances')
                return 0

        else:
            print('Please enter a meaningful 2D point group')
            return 0

    def rf_chir(self, dtype = 0):

        return self.sign*self.rf_grpchir(pgroup = 2, dtype = dtype)


    'Return oriented PF'
    def projform(self):

        a = sum(self.vec)

        return PF2([(self.r_02 - self.r_01)/a, (3*self.r_12)/a], self.sign)


    'Return Coform'
    def coform2d(self):

        return[i**2 for i in self.rightsign()]

    'Return Voform'
    def voform2d(self):

        c = self.coform2d()

        return[c[1]+c[2], c[0]+c[1], c[0]+c[2]]

    'Reconstitute lattice'
    def make2lat(self):

        l = self.coform2d()
        m = self.voform2d()
        cs = -l[0]/(np.sqrt(m[1])*np.sqrt(m[2]))
        alph = np.arccos(cs)

        v_1 = [np.sqrt(m[1]), 0]
        v_2 = [np.sqrt(m[2])*np.cos(alph), np.sqrt(m[2])*np.sin(alph)]

        return [v_1, v_2]

'----------------------------------------------------'
'PROJECTED FORM CLASS'
'----------------------------------------------------'
class PF2:

    def __init__(self, point, sign):
        self.qtpoint = point
        self.x = self.qtpoint[0]
        self.y = self.qtpoint[1]

        'Throws a warning if the projected form is not in QT'
        if self.x + self.y > 1 :
            print('Warning! Projected form is not in Quotient Triangle')

        'Chirality sign reverts to 0 if on boundary'
        if np.abs(1 - (self.x + self.y)) < tol  or (self.x < tol or self.y < tol):
            self.sign = 0
        else:
            self.sign = sign

    def qs_plot(self):
        if self.sign < tol:
            return [1-self.x, 1-self.y]
        else:
            return [self.x, self.y]

    'Calculates chirality based on infinity metric and position in quotient square'
    def pf_grpchir(self, dtype = 0, pgroup = 2):

        if dtype == 0:
            if pgroup == 2:
                return min([self.x, self.y, (1-self.x -self.y)/2])
            elif pgroup == 4:
                return self.x
            elif pgroup == 6:
                return (1-self.y)
            else:
                print('Please enter a meaningful point group!')
                return 0
        elif dtype == 2:
            if pgroup == 2:
                return min([self.x, self.y, (1-self.x -self.y)/2])
            elif pgroup == 4:
                return distgen(2, [self.x, self.y], [0,0])
            elif pgroup == 6:
                return distgen(2, [self.x, 1-self.y], [0,1])
            else:
                print('Please enter a meaningful point group!')
                return 0
        else:
            if pgroup == 2:
                print('Can only calculate D2 chirality for either L_2 or L_inf metric')
                return 0
            elif pgroup == 4:
                return distgen(dtype, [self.x, self.y], [0,0])
            elif pgroup == 6:
                return distgen(dtype, [self.x, 1-self.y], [0,1])
            else:
                print('Please enter a meaningful point group!')
                return 0

    def pf_chir(self, dtype = 0):

        return self.sign*self.pf_grpchir(pgroup = 2, dtype = dtype)

    def root_from_PF2(self):
        x = self.x
        y = self.y

        r_12 = y/3
        r_01 = (1-(r_12 + x))/2
        r_02 = (1 -r_12 + x)/2

        l = sorted([r_12, r_01, r_02])

        return RF2_signed(l, self.sign)

    def lattice_from_PF2(self):

        return self.root_from_PF2().make2lat()


'Calculate Chebyshev distance between two root forms. '
'Set orient = true to calculate oriented distance'
def rf2dist(rf_1, rf_2,  orient = True, dtype = 0):
    rfv_1 = rf_1.vec
    rfv_2 = rf_2.vec
    if not orient or ((rf_1.sign == rf_2.sign) or (rf_1.sign == 0 or rf_2.sign == 0)):
        return max(np.abs(rfv_1[0] - rfv_2[0]),np.abs(rfv_1[1] - rfv_2[1]),np.abs(rfv_1[2] - rfv_2[2]))
    else:
        if dtype == 0:
            d_0 = max(rfv_1[0] + rfv_2[0], np.abs(rfv_1[1]-rfv_2[1]), np.abs(rfv_1[2] - rfv_2[2]))
            d_1 = max(np.abs(rfv_1[2] - rfv_2[2]), minf(rfv_1[0], rfv_1[1], rfv_2[0], rfv_2[1]))
            d_2 = max(np.abs(rfv_1[0] - rfv_2[0]), minf(rfv_1[1], rfv_1[2], rfv_2[1], rfv_2[2]))
            return min(d_0, d_1, d_2)

        else:
            c1 = (-rfv_2[0], rfv_2[1], rfv_2[2])
            c2 = (rfv_2[2], rfv2_[0], rfv_2[1])
            c3 = (rfv_2[0], rfv_2[2], rfv_2[1])
            return min(distgen(2, rfv_1, c1), distgen(2, rfv_1, c2), distgen(2, rfv_1, c3))

'Calculate Chebyshev or L_2 distance between two projected forms. '
'Set orient = true to calculate oriented distance'
def pf2dist(pf_1, pf_2, orient = True, dtype =0):
    p_1 = pf_1.qtpoint
    p_2 = pf_2.qtpoint
    if not orient or (pf_1.sign == pf_2.sign or (pf_1.sign == 0 or pf_2.sign == 0)):
        return distgen (2, p_1, p_2)
    else:
        if dtype == 0:
            d_x = max([np.abs(p_2[0] - p_1[0]), p_2[1]+p_1[1]])
            d_y = max([p_2[0] + p_1[0], np.abs(p_2[1]-p_1[1])])
            d_xy = max([np.abs(p_2[0]-p_1[0]), 1-p_2[0]-p_2[1], np.abs(1-p_1[1]-p_2[0])])
            return min(d_x, d_y, d_xy)
        else:
            return min(distgen(2, p_1, [-p_2[0], p_2[1]]), distgen(2, p_1, [p_2[0], -p_2[1]]), distgen(2, p_1, [1-p_2[1], 1-p_2[0]]))
