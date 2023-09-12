# Rootform_2d
Calculation of and with root forms in 2 Dimensions

**Utility Functions**

*distgen(q, v_1, v_2)* is a generic Minkowkski distance calculator in R^2. It takes two inputs - a value of q >= 0 and two points in R^2 given as lists. If q=0 it will return the Chebyshev (L_inf) distance between points, otherwise it will return the L_q distance. 

*mincyc(l)* takes a list l of any length and returns the cyclic permutation that puts the smallest list member first

*roundlist(l,r)* takes a list of floats l and an integer value r. It returns l with r-rounded float values. 

*minf(a,b,c,d)* calculates the quantity max(|a-b|, |c-d|, |a+b-c+d|/2) required for the root form Chebyshev distance

*makelat(a,b,t)* will return the vectors of a Niggli Reduced 2D lattice given two lengths a, b and an angle t

*sb_sign(veclist)* will take a list of three vectors and, assuming that they are an obtuse superbase, calculate the sign of that superbase

*len_order(l)$ will take a list of vectors as arrays and order them by their length.  

*sb_sign(l)* will calculate the sign of an obtuse superbase from a list of three vectors (positive if superbase vectors are ordered anticlockwise, negative if superbase vectors are ordered clockwise)

**Lattice Generation and Manipulation**

*lat_from_obsb* will return the two shortest vectors of an input  obtuse superbase

*haar()* will generate a 'random' two dimensional lattice according to the Haar distribution. 

**Lat2D** - the 2D lattice class

The Lat2D class should be entered as two lists [[x_1, y_1], [x_2, y_2]]. The following methods are available:

*param_lat()* will return the length and angle parameters in the form [a,b,theta] for the lattice

*make_obsb()* will iterate the reduction step detailed in the paper until an obtuse superbase is reached and return the three vectors of that obtuse superbase as a list ordered by length

*make_cf* returns the **unsorted** coform for the lattice

*lattice_sign(tol = 10**-6)* returns the sign of the lattice as calculated from the coform - it reverses the sign of the superbase if a swap of two values is required to order the coform values list, or keeps that sign if ordering the list requires only permutation. It will the sign to 0 if any value in the coform is less than the input tolerance (10**-6 default). 

*make_rf(tol= 10**-16)* returns the Orientation-Aware root invariant object (see below) for the lattice- it will return 0 if the coform is achiral (to within a set floating point tolerance of $10^-16$), otherwise it will return the sign passed to it by the *lattice_sign* function. 

*make_pf(tol= 10**-16)* will return the Orientation-Aware projected invariant object (see below) for the lattice, again returning sign zero if it is within the set floating point tolerance level of bein

*make_qs(tol= 10**-16)* will return the relevant point in the quotient square for the signed lattice, calculated with the floating point tolerance as above

*build_lat(m, n)* will return a list of lattice points which are integer combinations av_1, bv_2, where a ranges between [-m,m] and b ranges between [-n, n]

*lattice_scale(self, scale_type = 'length')* will return a uniformly scaled lattice - it will either be such that v_1 = (1,0) if scale_type is set to 'length', or it will be a lattice of unit volume of scale_type is set to 'volume'. 

*mod_lat_points()* will return the three possible points for a lattice in the moduli space of lattices

*unit_cell_plot(file, title)* will save to the filename given in file a picture of the lattice unit cell with a title given by the 'title' string

*lat_plot(save_file, r = 2, type = 'Points', unit_cell = False, vectors = True, obsb = False)* will save to save_file an rxr grid plot of the lattice with thickened and arrowed lattice vectors in showing the unit cell. If unit_cell is set to True a filled unit cell will be displayed at the plots origin. If obsb is set to True the additional superbase vector will be shown. Type defaults to 'Points' which just shows lattice points. Set type to 'Voronoi' to view the Voronoi cells of the lattice, or to 'Delaunay' to view its Delaunay triangulation. 

*lat_to_SRI()* and $lat_to_SPI()* will return the spherical root and projected invariant of the lattice  respectively. 


**RI Class**

RI objects consist of a pair of inputs - an ordered list of three positive numbers and a sign (0, +1, -1). If an **achiral** root forms is entered, the sign will default to 0 regardless of the input. If the user inputs an unordered list there will be a warning message:

*rightsign()* will return the unordered root form - that is, it will swap the last two values of the root form if the sign is negative

*rf_grpchir(dtype = 0, pgroup = 2)* will, return the signed L_inf or L_2 based point group chirality of the root form. Id pgroup is set to 4 or 6, the D4 or D6 chirality will instead be calculated. If dtype is set to 2, the L_2 distance will be calculated rather than L_inf

*rf_chir(dtype = 0) will return the overall chirality for the root invariant. 

*projform()* will return the co-ordinates of the form in the quotient triangle, along with the sign of the root form, as a PF2 object. 


*make2lat()* will reconstitute the actual 2 dimensional lattice based on the input root form. Note that the output lattice will **not** be Niggli reduced - it will be the two vectors v_1, v_2 calculated as part of the obtuse superbase. 

**PI Class**

Projected Invariant (PI) objects should be entered  as a list [x, y] followed by a value +/-1 or 0 indicating the sign of the lattice from which the projected form is derived, and a tolerance set by defalut at 10**-6 - the sign will set to zero if either x or y value is within the tolerance level. The following methods are available:

*qs_plot()* will return the co-ordinates of the projected form in the quotient square. 

*pf_grpchir(dtype = 0, pgroup = 2)* returns the signed L_inf or L_2 chirality measure for the projected form depending on whether dtype is set to 0 or 2. Setting pgroup to 4 or 6 will calculate the D4 or D6 Chirality. Any other setting of dtype or pgroup will return an error message. 

*pf_chir(dtype = 0)* will return the signed D2 chirality. Defaults to use of Chebyshev distance - if dtype = 2 the Euclidean distane will be used instead. 

*root_from_PI (scale = 1)* will return a signed root form derived from an input projected form, which can be scaled by setting the scale factor

*lattice_from_PI (sc=1)* will return the 2D lattice reconstituted from an input projected form - setting the scale factor scales the underlying root form

*sphere_proj()* returns the latitude and longitude of the spherical projection. 

**SRI class**

Takes a spherical root invariant as a list of [latitude, longitude, size]

*xyz()* will return the co-ordinates of the point in R^3

*src_grp(grp)* will return the spherical root G-chiral distances. Set grp = 2, 4, 6 to compute D2, D4, or D6-chirality respectively

*src()* will give the spherical root chiral distance. 

*sri_to_lat()* will retrieve the lattice from the spherical root invariant

*sri_to_root()* will retrieve the root invariant from the spherical root invariant

**SPI class**
Takes a  spherical projected invariant as a list of [latitude, longitude]

*spi_to_pi(self)* will return the relevant projected invariant

*spi_to_lat(sc=1)* will reconstruct the vectors of a lattice of size given by sc from the SPI

*spc_grp(grp)* will return the spherical projected G-chiral distances. Set grp = 2, 4, 6 to compute D2, D4, or D6-chiral distance respectively

*spc()* will give the spherical projected chiral distance. 

**Distance Calculations**

*pf2dist(pf_1, pf_2, orient = True, dtype = 2)* will take two projected forms and calculate the L_2 oriented distance between them by default. It will calculate the unoriented distance if orient is set to False (or both projected forms share a sign or either are achiral), and the oriented or unoriented Chebyshev distance if dtype is set to 0

*rf2dist(rf_1, rf_2, orient = True, dtype = 2)* will take two root forms and calculate the L_2 oriented distance between them by default. It will calculate the unoriented distance if orient is set to False (or both projected forms share a sign or either are achiral), and the oriented or unoriented Chebyshev distance if dtype is set to 0

