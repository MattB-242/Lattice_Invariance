import numpy as np

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
        if l.index(i) = index:
            newlist.append(-eps)
        else:
            newlist.append(i+(2*eps))

    return newlist

'Cycle any list so that its least member is first'
def mincyc(l):

    top = l.index(min(l))

    return l[top:] + l[:top]


class RF2_signed:

    def __init__(self, vec, sign):
        self.unorient = vec
        self.sign = sign
        self.r_12 = vec[0]
        self.r_01 = vec[1]
        self.r_02 = vec[2]

    def rightsign(self):
        if self.sign = -1:
            return [self.unorient[0], self.unorient[2], self.unorient[1]]
        else:
            return self.unorient

    def projform(self):

        a = sum(vec)

        return [(self.sign*(self.r_01 - self.r_02))/a, self.r_12/a]


    def coform2d(self):

        return[i**2 for i in self.rightsign()]

    def voform2d(self):

        c = self.coform2d

        return[c[1]+c[2], c[0]+c[1], c[0]+c[2]]

    'Reconstitute lattice'
    def make2lat(self):

        l = self.coform2d()
        m = self.voform2d()
        alph = np.arccos(l[0]/(l[1]*l[2]))

        v_1 = [0,0, m[0]]
        v_2 = [np.cos(alph), self.sign*np.sin(alph)]

        return [v_1, v_2]

    def permlist(self, orient = True):
        if orient:
            return [self.unorient, [self.r_02, self.r_12, self.r_01], [self.r_01, self.r_02, self.r_12]]
        else:
            return [self.unorient, [self.r_02, self.r_12, self.r_01], [self.r_01, self.r_02, self.r_12],
                    self.r_01, self.r_12, self.r_02], [self.r_02, self.r_01, self.r_02], [self.r_12, self.r_02, self.r_01]]



'Calculate distance between two root forms'
def rf2dist(rf_1, rf_2, d = 2, orient = True):

    pass

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
        self.oblen = np.linalg.norm(ob)
        self.inner = np.dot(vec_1, vec_2)
        self.angle = np.rad2deg(np.arccos(self.inner/(self.xlen * self.ylen)))

    def make_CF(self):

        inners = [np.dot(self.x, self.y), np.dot(self.ob, self.x), np.dot(self.ob, self.y)]

        while min(inners) < 0:
            for i in inners:
                if i < 0:
                    inners = redstep(inners, inners.index(i))

        return mincyc(inners)

    def make_RF(self):

        c = self.make_CF()
        if c[2] > c[1]:
            return RF2_signed(np.sqrt(c_0), np.sqrt(c_2), np.sqrt(c_1)], -1)
        else:
            return RF2_signed(np.sqrt(c_0), np.sqrt(c_1), np.sqrt(c_2)], 1)

        
