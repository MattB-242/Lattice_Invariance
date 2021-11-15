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
        if l.index(i) == index:
            newlist.append(-eps)
        else:
            newlist.append(i+(2*eps))

    return newlist

'Cycle any list so that its least member is first'
def mincyc(l):

    top = l.index(min(l))

    return l[top:] + l[:top]

def roundlist(l):

    rl = []
    for i in l:
        rl.append(round(i,2))

    return rl

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
            return RF2_signed([np.sqrt(c_0), np.sqrt(c_2), np.sqrt(c_1)], -1)
        else:
            return RF2_signed([np.sqrt(c_0), np.sqrt(c_1), np.sqrt(c_2)], 1)


'----------------------------------------------------'
'2D ORIENTED ROOT FORM CLASS'
'----------------------------------------------------'

class RF2_signed:

    def __init__(self, vec, sign):
        self.unorient = vec
        self.r_12 = vec[0]
        self.r_01 = vec[1]
        self.r_02 = vec[2]
        if self.r_12 == 0 or (self.r_12 == self.r_01 or self.r_01 == self.r_02):
            self.sign = 1
        else:
            self.sign = sign

    'Create correctly ordered unoriented root form'
    def rightsign(self):
        if self.sign == -1:
            return [self.unorient[0], self.unorient[2], self.unorient[1]]
        else:
            return self.unorient

    'Calculate Chirality (L_inf ONLY)'
    def rf_chirality(self, dtype = 0):

        if dtype == 0:
            return min([self.r_12, (self.r_01 - self.r_12)/2, (self.r_02 - self.r_01)])
        else:
            print("Sorry, I don't know how to do that yet!")
            return 0

    'Find Nearest Chiral Lattice (L_inf ONLY)'
    def rf_nearest_chiral(self):

        pass

    'Return oriented PF'
    def projform(self):

        a = sum(vec)

        return [(self.sign*(self.r_01 - self.r_02))/a, self.r_12/a]


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

        v_1 = [m[0], 0]
        v_2 = [np.cos(alph), self.sign*np.sin(alph)]

        return [v_1, v_2]

    def permlist(self, orient = True):
        if orient:
            return [self.unorient, [self.r_02, self.r_12, self.r_01], [self.r_01, self.r_02, self.r_12]]
        else:
            return [self.unorient, [self.r_02, self.r_12, self.r_01], [self.r_01, self.r_02, self.r_12],
                    [self.r_01, self.r_12, self.r_02], [self.r_02, self.r_01, self.r_02], [self.r_12, self.r_02, self.r_01]]

'----------------------------------------------------'
'2D ORIENTED ROOT FORM CLASS'
'----------------------------------------------------'
class PF2:

    def __init__(self, point):
        self.x = point[0]
        self.y = point[1]
        if self.x < 0:
            self.sign = -1
        else:
            self.sign = 1

    def pf_chirality(self, dtype = 0):
        if dtype == 0:
            return min([self.x, self.y, (1-(2*self.x) - (3*self.y))/5])
        else:
            print("I don't know how to do that yet!")
            return 0


    def find_nearest_pf_point():

        pass

    def root_from_PF2(self):
        x = self.sign*self.x
        y = self.y

        r = (1 + 2*x - y)/2

        l = [y, 1 - (r + y), r]

        return RF2_signed(l, self.sign)

    def lattice_from_PF2(self):

        return self.root_from_PF2().make2lat()




'Calculate distance between two root forms'
def rf2dist(rf_1, rf_2, d = 0, orient = True):



    pass

'Calculate distance between two projected forms'
def rf2dist(rf_1, rf_2, d = 0, orient = True):

    pass


'PFs for Inverse Design'
p1 = [1/10, 1/10]
p2 = [1/6, 1/6]
p3 = [1/24, 1/6]
p4 = [1/4, 1/24]
p5 = [0,0]
p6 = [1/4, 1/6]
p7 = [0, 1/3]

allpfs = [p1, p2, p3, p4, p5, p6, p7]


for i in allpfs:
    pform = PF2(i)
    pchir = pform.pf_chirality()
    rform = pform.root_from_PF2()
    cform = rform.coform2d()
    vform = rform.voform2d()
    rchir = rform.rf_chirality()
    rlist = roundlist(rform.unorient)
    lat = [roundlist(i) for i in pform.lattice_from_PF2()]
    print('-------------------------------------------')
    print('for PF ' + repr(roundlist(i)))
    print('PC_Inf = ' + repr(round(pchir,2)))
    print('Root Form = ' + repr(rlist))
    print('Coform = ' + repr(roundlist(cform)))
    print('Voform = ' + repr(roundlist(vform)))
    print('RC_Inf =' + repr(round(rchir,2)))
    print('Actual Lattice = ' + repr(lat))
    print('-------------------------------------------')
