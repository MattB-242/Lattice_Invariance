import numpy as np
import itertools as it

def main():

    testlist_a = [4,2,9]
    testlist_b = [1,4,8]

    testsquare = Lattice_2d([[1,0], [-0.5, 0.85]])
    testhex = Lattice_2d([[1,0], [-0.5, 0.75]])

    testcube = Lattice_3d([[1,0,0], [0,1,0], [0,0,1]])
    testother = Lattice_3d([[1,0,0], [0.1,1,0], [0.1,0.2,1]])

    squobase = Obsuper_2d(testsquare.msb())
    hxobase = Obsuper_2d(testhex.msb())


    print('Distance between hexagonal and square base')
    print(minperdist(squobase.sellings_t,hxobase.sellings_t, 1))

    cubase = Obsuper_3d(testcube.msb3d())
    obase = Obsuper_3d(testother.msb3d())


    print('Distance between Cube and Test base')
    print(minperdist(cubase.sellings_h,obase.sellings_h, 1))

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


def minperdist (a,b, l=1, cyc = False):

        if len(a) != len(b):
            print('Both lists must be the same length')
            return 0

        elif cyc:
            perlist = [[a[j-i] for i in range (len(a))] for j in range (len(a))]
        
        else:
            perlist = list(it.permutations(a))

            allperdist = []
            minlists = []
            for i in perlist:
                allperdist.append(dst(i,b,l))
                minlists.append(i)

        print('minimal L_'+repr(l) + ' distance.')
        best = min(allperdist)
        return best

def angcalc_2d(a=1, b=1, alpha):
    
    return [[a,0], [b*np.cos(alpha), b*np.sin(alpha)]
    

class Lattice_2d:

    def __init__(self,lattice):
        self.lattice = lattice
        self.x = np.array(self.lattice[0])
        self.y = np.array(self.lattice[1])
        self.inner = np.dot(self.x, self.y)
        self.angle = np.arccos(self.inner)
        self.mat = np.transpose(np.array(self.lattice))
        self.vol = np.linalg.det(self.mat)

    def msb(self):

        return[-np.add(self.x, self.y), self.x, self.y]

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

class Lattice_3d:

    def __init__(self,lattice):
        self.lattice = lattice
        self.x = np.array(lattice[0])
        self.y = np.array(lattice[1])
        self.z = np.array(lattice[2])
        #self.mat = np.array([[[lattice[0][0], lattice[1][0]], [[lattice[0][1],lattice[1][1]]])
        #self.vol = np.linalg.det(self.mat)

    def msb3d(self):

        return[-np.add(self.x, self.y, self.z), self.x, self.y, self.z]

class Obsuper_3d:

    def __init__(self,basis):
        self.basis = basis
        self.vo = np.array(basis[0])
        self.vl = np.array(basis[1])
        self.vt = np.array(basis[2])
        self.vh = np.array(basis[3])
        self.p_ol = np.dot(self.vo, self.vl)
        self.p_ot = np.dot(self.vo, self.vt)
        self.p_oh = np.dot(self.vo, self.vh)
        self.p_lt = np.dot(self.vl, self.vt)
        self.p_lh = np.dot(self.vl, self.vh)
        self.p_th = np.dot(self.vt, self.vh)
        self.sellings_h = [0, -self.p_ol, -self.p_ot, - self.p_oh, -self.p_lt, -self.p_lh, -self.p_th]

if __name__ == "__main__":
    main()
