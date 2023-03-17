import numpy as np
import pandas as pd
import csv
import random
import math
import matplotlib
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d
from scipy.spatial import Delaunay, delaunay_plot_2d
from sklearn.metrics.pairwise import haversine_distances as hv

matplotlib.rcParams['text.usetex'] = True

'-----------------------------------------------------'
'UTILITY FUNCTIONS'
'-----------------------------------------------------'
'General L_q distance calculator for vectors in R_2'
def distgen(q, v_1, v_2):

    if q == 0:
        return max(np.abs(v_1[0] - v_2[0]), np.abs(v_1[1] - v_2[1]))
    else:
        return ((np.abs(v_1[0] - v_2[0])**q) + (np.abs(v_1[1] - v_2[1])**q))**(1/q)

'Cycle any list so that its least member is first'
def mincyc(l):

    top = l.index(min(l))

    return l[top:] + l[:top]

'Round a list of floats to some defined level (default = 2dp)'
def roundlist(l, r=2):

    rl = []
    for i in l:
        rl.append(round(i,r))

    return rl

'Calculate M_inf quantity (see maths paper, Lemma C.2)'
def minf(a,b,c,d):
    return max([np.abs(a-b), np.abs(c-d), (np.abs(a+b-c+d))/2])

'Make a pair of 2D Lattice  vectors from length and angle parameters (in degrees)'
def makelat (a, b, t):
    trad = np.deg2rad(t)

    return Lat2d([a, 0], [b*np.cos(trad), b*np.sin(trad)])


'Native Python argsort function:'
'Given a list of quantities, return a list consisting of the index each quantity'
'would be at in the sorted version of the list'
def index_sorted(l):

    scomp = sorted(l.copy())
    indl = []
    countlist = []
    for i in l:
        shift = countlist.count(i)
        indl.append(scomp.index(i) + shift)
        countlist.append(i)

    return indl

'Order a list of vectors by length'
def len_order(l):
    poslist = list(np.argsort([np.linalg.norm(i) for i in l]))
    ordlist = []
    for i in range(0,max(poslist)+1):
        ordlist.append(l[poslist.index(i)])

    return ordlist


'Given an obtuse superbase as a list of vectors, calculate the lattice sign from determinant'
def sb_sign(l):
    latm = np.array(len_order(l)[:2])

    latang = np.abs(np.dot(latm[0], latm[1]))
    latlendiff = np.abs(np.linalg.norm(latm[1]) - np.linalg.norm(latm[0]))


    if latang == 0 or latlendiff == 0:
        return 0
    elif np.linalg.det(latm) < 0:
        return -1
    else:
        return 1

def lat_from_obsb(obsb):

    vecs = len_order(obsb[0])[:2]

    return(Lat2d(vecs[0], vecs[1]))

'Calculate a projected form from a longitude and latitiude'
def globe_project(mu, phi):
    if -157.5 <= mu < 180:
        mu_p = np.deg2rad(mu + 157.5)
    else:
        mu_p = np.deg2rad(517.5 + mu)
    mu_p = np.deg2rad(mu + 157.5)
    phi_p = (90 - np.abs(phi))/90
    r = 1 - (1/np.sqrt(2))

    if (np.rad2deg(mu_p) <= 112.5) or (np.rad2deg(mu_p) > 337.5):
        x =  (phi_p*(1-(2*r))*np.cos(mu_p))/(np.sin(mu_p) + np.cos(mu_p))
        y = (phi_p*(1-(2*r))*np.sin(mu_p))/(np.sin(mu_p) + np.cos(mu_p))

        if phi > 0:
            return PF2([x+r,y+r], 1)
        elif phi < 0:
            return PF2([x+r,y+r], -1)
        else:
            return PF2([x+r,y+r],0)

    elif 112.5 < np.rad2deg(mu_p) < 225:
        x = -r*phi_p
        y = -r*(phi_p * (np.sin(mu_p)/np.cos(mu_p)))

        if phi > 0:
            return PF2([x+r,y+r], 1)
        if phi < 0:
            return PF2([x+r,y+r], -1)
        else:
            return PF2([x+r,y+r],1)

    elif 225 <= np.rad2deg(mu_p) <= 337.5:
        y = -r*phi_p
        x = -r*(phi_p * (np.cos(mu_p)/np.sin(mu_p)))

        if phi > 0:
            return PF2([x+r,y+r], 1)
        elif phi < 0:
            return PF2([x+r,y+r], -1)
        else:
            return PF2([x+r,y+r],1)

'Generate a Haar Random Lattice'
def haar():
    x_1 = random.uniform(-0.5,0.5)
    x_2 = random.uniform(-0.5,0.5)
    x_3 = random.uniform(-0.5,0.5)

    u = np.sin((x_1*np.pi/3))
    v = (np.cos((x_1*np.pi/3)))/(0.5-x_2)
    w = np.pi*x_3

    a = np.array([[1,u], [0,1]])
    b = np.array([[np.sqrt(v), 0],[0, 1/np.sqrt(v)]])
    ab = np.matmul(a,b)

    c = np.array([[np.cos(w), -np.sin(w)],[np.sin(w), np.cos(w)]])

    fin = np.matmul(ab,c)

    return Lat2d(fin[0], fin[1])

'Add a row of lattice invariant data to a file using either the PI, spherical co-ordinates or lattice parameters'
def row_input(fle, given, pair, name = 'None'):

    df = pd.read_csv(fle, sep = ',')

    if given == 'spherical':
        pform = globe_project(pair[0], pair[1])
        lat = pform.lattice_from_PF2()
        param_lat = lat.param_lat()
        ang = param_lat[2]
        len = param_lat[1]/param_lat[0]


        insert_row = {'name':name, 'long':pair[0], 'lat':pair[1], 'x':pform.x,'y':pform.y, 'v2':len, 'theta':ang}

    elif given == 'projected':
        if pair[0] + pair[1] > 1:
            pform = PF2(pair, -1)
        else:
            pform = PF2(pair, 1)

        proj = pform.sphere_proj()

        lat = pform.lattice_from_PF2()
        param_lat = lat.param_lat()
        ang = param_lat[2]
        len = param_lat[1]/param_lat[0]

        insert_row = {'name':name, 'long':proj.lon, 'lat':proj.lat, 'x':pair[0], 'y':pair[1], 'v2':len, 'theta':ang}

    elif given == 'parameter':

        lat = makelat(1,pair[0], pair[1])
        pform = lat.make_pf()
        proj = pform.sphere_proj()

        insert_row = {'name':name, 'long':proj.lon, 'lat':proj.lat, 'x':pform.x, 'y':pform.y, 'v2':pair[0], 'theta':pair[1]}

    else:
        print('Input data type not recognised')

        insert_row = {'name': 0, 'long': 0, 'lat': 0, 'x': 0, 'y': 0, 'v2': 0,'theta': 0}

    print('Inserting row')
    print(insert_row)


    df = pd.concat([df, pd.DataFrame([insert_row])], ignore_index = True)
    if 'Unnamed: 0' in df.columns.values.tolist():
        df.drop('Unnamed: 0', axis=1, inplace=True)



'----------------------------------------------------'
'LATTICE CLASS IN 2D WITH ROOT FORM GENERATOR'
'----------------------------------------------------'

'Basic Lattice Class - takes two input vectors in R2'
class Lat2d:

    def __init__(self, vec_1, vec_2):
        self.veclist = [vec_1, vec_2]
        self.x = np.array(vec_1)
        self.y = np.array(vec_2)
        self.ob = -(self.x + self.y)
        self.xlen = np.linalg.norm(vec_1)
        self.ylen = np.linalg.norm(vec_2)
        self.oblen = np.linalg.norm(self.ob)
        self.inner = np.dot(vec_1, vec_2)
        self.angle = np.rad2deg(np.arccos(self.inner/(self.xlen * self.ylen)))
        self.volume = np.linalg.det(np.array([len_order([vec_1, vec_2])]))
        self.scale = np.abs(self.volume)

    'Return a lattice in terms of its length and angle parameters'
    def param_lat(self):

        return([self.xlen, self.ylen, self.angle])

    'Reduce any input lattice to its obtuse superbase'
    def make_obsb(self, step = False):

        v_1 = self.x
        v_2 = self.y
        v_0 = self.ob


        inners = [-np.dot(v_1, v_2), -np.dot(v_0, v_1), -np.dot(v_0, v_2)]

        stepcount = 0

        if min(inners) >= 0:

            return ([[v_0, v_1, v_2], stepcount])

        else:
            while min(inners) < 0:
                if step:
                    print('Step ' + repr(stepcount) +':')
                    print('Superbase vectors: ' + repr([[v_0, v_1,v_2]]))
                    print('Inner products: ' + repr(inners))
                    input('Press Enter to Run a Reduction Step: ')
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

                stepcount +=1

            return ([[v_0, v_1, v_2], stepcount])

    'Return the two shortest vectors of the obtuse superbase of the lattice'


    'Generate the Coform of the obtuse superbase of a lattice as a list'
    def make_cf(self):

        obsb = self.make_obsb()[0]

        return [-np.dot(obsb[1], obsb[2]), -np.dot(obsb[0], obsb[1]), -np.dot(obsb[0], obsb[2])]


    'Calculate the sign of a lattice'
    def lattice_sign(self, tol=10**-6):

        return sb_sign(self.make_obsb()[0])

    'Calculate the root form of a lattice - returns a ROOT FORM object'
    def make_rf(self, tol=10**-16):

        cf = self.make_cf()
        #print('we start with the coform ' + repr(cf))
        sgn = self.lattice_sign()
        #print('The lattice has sign ' + repr(sgn))

        rf = sorted([np.sqrt(i) for i in cf])
        #print('now the root form is ' + repr(rf))

        if np.abs(rf[0])<tol or (np.abs(rf[0] - rf[1]) < tol or np.abs(rf[1] - rf[2])<tol):
                return RF2_signed(rf, 0)
        else:
                return RF2_signed(rf, sgn)

    'Calculate the projected form of a lattice'
    def make_pf(self, tol = 10**-16):

        return self.make_rf().projform()

    'Calculate the position of a lattice in the quotient square'
    def make_qs(self, tol = 10**-6):

        return self.make_rf().projform().qs_plot()

    'Given a lattice basis v_1, v_2, return all lattice points of the form av_1 + bv_2'
    'for a in the range [-m,m] and b in the range [-n,n]'
    def build_lat(self, m, n):

        '''
        if self.x[1] != 0:
            print('Warning! The input lattice is not rotated to the x-axis.')
        '''

        layer_list = []

        for i in range(-n,n+1):
            for j in range (-m,m+1):
                layer_list.append(i*self.x + j*self.y)

        return layer_list

    'Pin lattice to the x-axis'
    def axis_rotate(self):

        sheep_dip = self.make_pf().lattice_from_PF2()

        return Lat2d(sheep_dip.x, sheep_dip.y)


    'Rotate and scale any 2D lattice to either det 1 or |x| = 1'
    def lattice_scale(self, scale_type = 'length'):
        k = True

        if self.x[1] == 0: #Rotate lattice to x-axis if not already there
            rotlat = self
        else:
            rotlat = self.axis_rotate()

        if scale_type == 'length': #Generate scaling factor
            sc = 1/rotlat.x[0]
            return Lat2d([self.x[0] * sc, 0], [self.y[0] * sc, self.y[1] * sc])
        elif scale_type == 'volume':
            sc = 1/(np.det(np.array([self.x, self.y])))**2
            return Lat2d([self.x[0] * sc, 0], [self.y[0] * sc, self.y[1] * sc])
        else:
            print('Error:scaling type not recognised')
            print('Original lattice returned unchanged')
            return self

    'Return a list of moduli space points representing the input lattice'
    def mod_lat_points(self):

        def refl(v):
            return [-1-v[0], v[1]]

        def inv(v):
            c = v[0]**2 + v[1]**2

            return [v[0]/c, v[1]/c]

        def shift(v):
            return[v[0]+1, v[1]]

        inlat = self.lattice_scale(scale_type = 'length')
        c = 1/inlat.x[0]
        v_2 = [c*inlat.y[0], c*inlat.y[1]]


        pts =  [refl(v_2), inv(v_2), inv(refl(inv(v_2)))]

        modpts = []

        for i in pts:
            if i[0] < -0.5:
                modpts.append(shift(i))
            else:
                modpts.append(i)

        return modpts

    def unit_cell_plot(self, file, title):

        unit_cell_x = [0, self.x[0], self.x[0] + self.y[0], self.y[0]]
        unit_cell_y = [0, self.x[1], self.x[1] + self.y[1], self.y[1]]

        plt.scatter(unit_cell_x, unit_cell_y, zorder = 4)
        plt.xticks(color = 'white')
        plt.yticks(color = 'white')
        plt.tick_params(bottom = False)
        plt.tick_params(left = False)
        ax = plt.gca()
        ax.set_aspect('equal', adjustable='box')
        plt.arrow(0,0,self.x[0], self.x[1], width = 0.01, color = 'black', length_includes_head = True, zorder = 3)
        plt.arrow(0,0,self.y[0], self.y[1], width = 0.01, color = 'black', length_includes_head = True, zorder = 3)
        plt.plot(unit_cell_x, unit_cell_y,linewidth = 2, color="green", zorder = 2)
        plt.fill(unit_cell_x, unit_cell_y, color = 'yellow', zorder = 1)
        plt.title(title, fontsize = 25)
        plt.box(False)
        plt.savefig(file)
        #plt.show()
        plt.cla()

    'Create an r x r grid of the lattice with or without vectors and a filled unit cell.'
    def lat_plot(self, save_file, r = 2, type = 'Points', unit_cell = False, vectors = True):
        latpoints = []
        unit_cell_x = [0, self.x[0], self.x[0] + self.y[0], self.y[0]]
        unit_cell_y = [0, self.x[1], self.x[1] + self.y[1], self.y[1]]
        #print('unit cell coordinates {}.'.format([unit_cell_x,unit_cell_y]))

        #print('making base layer')
        basepoints = []
        for i in range(-r,r+1):
            basepoints.append(i*np.array(self.x))
            latpoints += basepoints
            #print('making upper and lower layers')
            for j in range (-r, r+1):
                #print(j)
                if j!= 0:
                    layerpoints = []
                    #print('making a new layer')
                    for k in basepoints:
                        #print('adding to {}'.format(k))
                        layerpoints.append(k + j*np.array(self.y))
                        latpoints += layerpoints


        if type == 'Points':
            x = [i[0] for i in latpoints]
            y = [i[1] for i in latpoints]

            plt.scatter(x,y, zorder = 3)
            plt.xticks(color = 'white')
            plt.yticks(color = 'white')
            plt.tick_params(bottom = False)
            plt.tick_params(left = False)
            ax = plt.gca()
            ax.set_aspect('equal', adjustable='box')
            if vectors:
                plt.arrow(0,0,self.x[0], self.x[1], width = 0.01, color = 'black', length_includes_head = True, zorder = 2)
                plt.arrow(0,0,self.y[0], self.y[1], width = 0.01, color = 'black', length_includes_head = True, zorder = 2)
                plt.plot(unit_cell_x, unit_cell_y,linewidth = 2, color="green", zorder = 2)
            if unit_cell:
                plt.fill(unit_cell_x, unit_cell_y, color = 'yellow', zorder = 1)
            plt.box(False)
            plt.savefig(save_file)
            plt.show()
            plt.cla()

        elif type == 'Voronoi':
            vor = Voronoi(latpoints)
            vorplot = voronoi_plot_2d(vor, zorder = 3)
            ax = vorplot.gca()
            ax.set_aspect('equal', adjustable='box')
            ax.xaxis.set_ticklabels([])
            ax.xaxis.set_ticks([])
            ax.yaxis.set_ticklabels([])
            ax.yaxis.set_ticks([])
            plt.tick_params(bottom = False)
            plt.tick_params(left = False)
            if vectors:
                plt.arrow(0,0,self.x[0], self.x[1], width = 0.005, color = 'red', zorder = 2)
                plt.arrow(0,0,self.y[0], self.y[1], width = 0.005, color = 'red', zorder = 2)
            if unit_cell:
                plt.fill(unit_cell_x, unit_cell_y, color = 'yellow', zorder = 1)
            plt.box(False)
            plt.savefig(save_file)
            plt.show()
            plt.cla()

        elif type == 'Delaunay':
            dl = Delaunay(latpoints)
            dlplot = delaunay_plot_2d(dl)
            ax = dlplot.gca()
            ax.set_aspect('equal', adjustable='box')
            ax.xaxis.set_ticklabels([])
            ax.xaxis.set_ticks([])
            ax.yaxis.set_ticklabels([])
            ax.yaxis.set_ticks([])
            plt.tick_params(bottom = False)
            plt.tick_params(left = False)
            plt.savefig(save_file)
            plt.box(False)
            plt.show()
            plt.cla()

        else:
            print('Type not recognised!')

    'Retrieve Sphere Invariant from lattice'
    def lat_to_SRI(self):

        lat_h = sum(self.make_rf().vec)
        sphere_h = self.make_pf().sphere_proj()


        return SRI(sphere_h.lat, sphere_h.lon, lat_h)

    'Plot Lattice SLP'
    def lat_to_SLP(self):

        return(self.lat_to_SRI().spi_to_R2())


'----------------------------------------------------'
'2D ORIENTED ROOT FORM CLASS'
'----------------------------------------------------'

'Basic Root Form Class - takes a list of positive numbers and a sign'
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

    'Calculate Positive G-Chirality for groups D2, D4, D6 (L_inf and L_2)'
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
                return min([self.r_12,  (self.r_02-self.r_01)/2])
            elif dtype == 2:
                return np.sqrt((self.r_12)**2 + 0.25*(self.r_02 - self.r_02)**2)

        elif pgroup == 6:
            if dtype == 0:
                return (self.r_02 - self.r_12)/2
            elif dtype == 2:
                return np.sqrt(2/3 *(self.r_12**2 + self.r_01**2 + self.r_02**2 -(self.r_01*self.r_12) - (self.r_12*self.r_02) - (self.r_01*self.r_02)))


        else:
            print('Please enter a meaningful 2D point group')
            return 0

    'Calculate the signed D2 chirality (i.e. overall chirality)'
    def rf_chir(self, dtype = 0):

        return self.sign*min([self.rf_grpchir(pgroup = 2, dtype = dtype), self.rf_grpchir(pgroup = 4, dtype = dtype), self.rf_grpchir(pgroup = 6, dtype = dtype)])


    'Return projected form'
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

        return Lat2d(v_1, v_2)


'----------------------------------------------------'
'PROJECTED FORM CLASS'
'----------------------------------------------------'

'Projected Form Basic Class - takes point (x,y) as list and a sign'
class PF2:

    def __init__(self, point, sign, tol = 5*10**-8):
        self.qtpoint = point
        self.x = self.qtpoint[0]
        self.y = self.qtpoint[1]
        self.out_of_bounds = False

        'Throws a warning if the projected form is not in QT'
        if self.x + self.y > 1 + tol :
            self.out_of_bounds = True
            #print('Warning! Projected form is not in Quotient Triangle')

        'Chirality sign reverts to 0 if on boundary'
        if np.abs(1 - (self.x + self.y)) < tol  or (self.x < tol or self.y < tol):
            self.sign = 0
        else:
            self.sign = sign

    'Plots co-ordinates in the quotient square'
    def qs_plot(self):
        if self.sign == 0:
            return [self.x, self.y]
        elif self.sign < 0:
            return [1-self.y, 1-self.x]
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
                return distgen(2, [self.x, self.y], [0,1])
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
                return distgen(dtype, [self.x, self.y], [0,1])
            else:
                print('Please enter a meaningful point group!')
                return 0

    def pf_chir(self, dtype = 0):

        return self.sign*min([self.pf_grpchir(pgroup = 2, dtype = dtype), self.pf_grpchir(pgroup = 4, dtype = dtype), self.pf_grpchir(pgroup = 6, dtype = dtype)])

    'Returns root form from projected form at a given scale '
    def root_from_PF2(self, scale = 1):
        x = self.x
        y = self.y

        r_12 = scale * (y/3)
        r_01 = scale *((1-(r_12 + x))/2)
        r_02 = scale *((1 -r_12 + x)/2)

        l = sorted([r_12, r_01, r_02])

        return RF2_signed(l, self.sign)

    def lattice_from_PF2(self, sc = 1):

        return self.root_from_PF2(scale = sc).make2lat()

    'Plots spherical root invariant based on projected form co-ordinates'
    'Set scale to the height of the lattice if this is available, otherwise defaults to 1'
    def sphere_proj(self, sc = 1):
        t = 1-(1/np.sqrt(2))

        if self.x == t:
            phi = self.sign*90
            mu = np.nan
        else:
            psi = np.rad2deg(math.atan((self.y-t)/(self.x-t)))


            'calculate longitude'
            if self.x < t:
                mu = psi + 22.5
            elif self.x > t:
                if psi >= -22.5:
                    mu = psi-157.5
                else:
                    mu = psi+202.5

            'calculate latitude'
            if -180 < mu < -45:
                phi = self.sign*((1-self.x - self.y)/(np.sqrt(2)-1))*90
                if np.abs(phi) > 90:
                    print('hypotenuse problem!')

            elif -45 <= mu < 67.5:
                phi = self.sign*((np.sqrt(2)*self.x)/(np.sqrt(2)-1))*90
                if np.abs(phi) > 90:
                    print('vertical problem!')

            else:
                phi = self.sign*((np.sqrt(2)*self.y)/(np.sqrt(2)-1))*90
                if np.abs(phi) > 90:
                    print('horizontal problem!')

        return SRI(phi, mu, sc)

'----------------------------------------------------'
'LATTICE DISTANCE CALCULATIONS'
'----------------------------------------------------'

'Calculate Chebyshev distance between two lattices. '
'Set orient = true to calculate oriented distance'
def rf2dist(l_1, l_2,  orient = True, dtype = 0):
    rfv_1 = l_1.make_rf().vec
    rfv_2 = l_2.make_rf().vec
    if not orient or ((rfv_1.sign == rfv_2.sign) or (rfv_1.sign == 0 or rfv_2.sign == 0)):
        return max(np.abs(rfv_1[0] - rfv_2[0]),np.abs(rfv_1[1] - rfv_2[1]),np.abs(rfv_1[2] - rfv_2[2]))
    else:
        if dtype == 0:
            d_0 = max(rfv_1[0] + rfv_2[0], np.abs(rfv_1[1]-rfv_2[1]), np.abs(rfv_1[2] - rfv_2[2]))
            d_1 = max(np.abs(rfv_1[2] - rfv_2[2]), minf(rfv_1[0], rfv_1[1], rfv_2[0], rfv_2[1]))
            d_2 = max(np.abs(rfv_1[0] - rfv_2[0]), minf(rfv_1[1], rfv_1[2], rfv_2[1], rfv_2[2]))
            return min(d_0, d_1, d_2)

        else:
            c1 = (-rfv_2[0], rfv_2[1], rfv_2[2])
            c2 = (rfv_2[2], rfv_2[0], rfv_2[1])
            c3 = (rfv_2[0], rfv_2[2], rfv_2[1])
            return min(distgen(2, rfv_1, c1), distgen(2, rfv_1, c2), distgen(2, rfv_1, c3))

'Calculate Chebyshev or L_2 distance between two lattices. '
'Set orient = true to calculate oriented distance'
def pf2dist(l_1, l_2, orient = True, dtype =0):
    p_1 = l_1.make_pf().qtpoint
    p_2 = l_2.make_pf().qtpoint
    if not orient or (p_1.sign == p_2.sign or (p_1.sign == 0 or p_2.sign == 0)):
        return distgen (2, p_1, p_2)
    else:
        if dtype == 0:
            d_x = max([np.abs(p_2[0] - p_1[0]), p_2[1]+p_1[1]])
            d_y = max([p_2[0] + p_1[0], np.abs(p_2[1]-p_1[1])])
            d_xy = max([np.abs(p_2[0]-p_1[0]), 1-p_2[0]-p_2[1], np.abs(1-p_1[1]-p_2[0])])
            return min(d_x, d_y, d_xy)
        else:
            return min(distgen(2, p_1, [-p_2[0], p_2[1]]), distgen(2, p_1, [p_2[0], -p_2[1]]), distgen(2, p_1, [1-p_2[1], 1-p_2[0]]))

'Calculate Haversine Distance Between Two 2D Lattice Objects'
def lat_haver(l,m,r=1):

    sp_1 = l.lat_to_SRI()
    sp_2 = m.lat_to_SRI()

    return r*hv([sp_1,sp_2])[0][1]

def chord_3d(l,m):

    h_1 = sum(l[0].make_rf().vec)
    h_2 = sum(l[1].make_rf().vec)

    sp_1 = [np.deg2rad(i) for i in l[0].make_pf().sphere_proj().reverse()]
    sp_2 = [np.deg2rad(i) for i in l[1].make_pf().sphere_proj().reverse()]

'----------------------------------------------------'
'SPHERICAL ROOT INVARIANT CLASS'
'----------------------------------------------------'

'SRI - takes latitude, longitude and height'

class SRI:

    def __init__(self, lat, lon, h):
        self.lat = lat
        self.lon = lon
        self.rlat = np.deg2rad(lat)
        self.rlon = np.deg2rad(lon)
        self.height = h

    'Calculate xyz co-ordinates'
    def xyz(self):

        return np.array([self.height*np.cos(self.lat)*np.cos(self.lon), self.height*np.cos(self.lat)*np.sin(self.lon), self.height*np.sin(self.lat)])

    'Calculate Spherical Root Chirality'
    def src_grp(self,grp):

        if grp == 2:
            return self.height*2*np.abs(np.sin(np.deg2rad(self.lat/2)))

        elif grp == 4:
            return self.height * np.sqrt((np.cos(self.rlat)*np.cos(self.rlon) - np.cos(np.deg2rad(67.5)))**2 +
                                         (np.cos(self.rlat)*np.sin(self.rlon) - np.sin(np.deg2rad(67.5)))**2 +
                                         (np.sin(self.rlat))**2)

        elif grp == 6:
            return self.height * np.sqrt((np.cos(self.rlat)*np.cos(self.rlon) - np.cos(np.deg2rad(45)))**2 +
                                         (np.cos(self.rlat)*np.sin(self.rlon) + np.sin(np.deg2rad(45)))**2 +
                                         (np.sin(self.rlat))**2)

    'Calculate Spherical Root Chirality'
    def src(self):

        return self.src_grp(2)

    'Calculate Spherical Projected Group Chirality'
    def spc_grp(self, grp):

        if grp == 2:

            return np.abs(self.lat)

        elif grp == 4:

            return np.rad2deg(([[np.deg2rad(self.lat), np.deg2rad(self.lon)], [0, np.deg2rad(67.5)]])[0][1])

        elif grp == 6:

            return np.rad2deg(hv([[np.deg2rad(self.lat), np.deg2rad(self.lon)], [0, np.deg2rad(-45)]])[0][1])

        else:
            print('Error: Group should be 2, 4 or 6!')

    'Calculate Spherical Projected Chirality'
    def spc(self):

        return self.spc_grp(2)

    'Retrieve a projected invariant from an SRI'
    def spi_to_proj(self):

        return(globe_project(self.lon, self.lat))

    'Retrieve a lattice from an SRI'
    def spi_to_lat(self, sc = 1):

        return(globe_project(self.lon, self.lat).lattice_from_PF2(sc = self.height))

    'Return a spherical projection'
    def spi_to_R2(self):

        return [np.sin(self.rlon)*np.cos(self.rlat)/(1+(np.cos(self.rlon)*np.cos(self.rlat))),
                np.sin(self.rlat)/(1+(np.cos(self.rlon)*np.cos(self.rlat)))]
