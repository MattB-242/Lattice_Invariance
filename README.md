# Rootform_2d
Calculation of and with root forms in 2 Dimensions

**Utility Functions**

*distgen(q, v_1, v_2)* is a generic Minkowkski distance calculator in R^2. It takes two inputs - a value of q >= 0 and two points in R^2 given as lists. If q=0 it will return the Chebyshev (L_inf) distance between points, otherwise it will return the L_q distance. 

*mincyc(l)* takes a list l of any length and returns the cyclic permutation that puts the smallest list member first

*roundlist(l,r) * takes a list of floats l and an integer value r. It returns l with r-rounded float values. 

**Working with Root Forms**

Root Form objects consist of a pair of inputs - an ordered list of three positive numbers and a sign (0, +1, -1). If an **achiral** root forms is entered, the sign will default to 0 regardless of the input:

*rightsign* will return the unordered root form - that is, it will swap the last two values of the root form if the sign is negative

*rf_chir(dtype = 0, pgroup = 2)* will, return the signed L_inf or L_2 based point group chirality of the root form. Id pgroup is set to 4 or 6, the D4 or D6 chirality will instead be calculated. If dtype is set to 2, the L_2 distance will be calculated rather than L_inf

*projform* will return the co-ordinates of the form in the quotient triangle, along with the sign of the root form, as a PF2 object. 

*coform2d* and *voform2d* will return the coform  and voform for the obtuse superbase represented by the root form. Vectors are in anticlockwise order as stated in the paper. 

*make2lat* will reconstitute the actual 2 dimensional lattice based on the input root form. 

**Working with Projected Forms**

Projected Form (PF2) objects should be entered  as a list [x, y] followed by a value +/-1 or 0 indicating the sign of the lattice from which the projected form is derived. The following methods are available:

*qs_plot* will return the co-ordinates of the projected form in the quotient square. 

*pf_chir(dtype = 0, pgroup = 2)* returns the signed L_inf or L_2 chirality measure for the projected form depending on whether dtype is set to 0 or 2. Setting pgroup to 4 or 6 will calculate the D4 or D6 Chirality. Any other setting of dtype or pgroup will return an error message. 

*root_from_pf2* will return the signed root form derived from an input projected form

*lattice_from_pf2* will return the 2D lattice reconstituted from an input projected form 

**Working with Lattices**

The Lat2D class should be entered as two lists [[x_1, y_1], [x_2, y_2]]. The following methods are available:

*make_obsb* will iterate the reduction step detailed in the paper until an obtuse superbase is reached

*sb_sign* will calculate the sign of the superbase (positive if superbase vectors are ordered anticlockwise, negative if superbase vectors are ordered clockwise)

*make_cf* returns the **unsorted** coform for the lattice

*lattice_sign* returns the sign of the lattice as calculated from the coform - it reverses the sign of the superbase if a swap of two values is required to order the coform values list, or keeps that sign if ordering the list requires only permutation. 

*make_rf* returns the signed root form object - it will return 0 if the coform is achiral (to within a set tolerance of $10^-8$), otherwise it will return the sign passed to it by the *lattice_sign* function. 

**Distance Calculations**

*pf2dist(pf_1, pf_2, orient = True, dtype = 2)* will take two projected forms and calculate the L_2 oriented distance between them by default. It will calculate the unoriented distance if orient is set to False (or both projected forms share a sign or either are achiral), and the oriented or unoriented Chebyshev distance if dtype is set to 0

*rf2dist(rf_1, rf_2, orient = True, dtype = 2)* will take two root forms and calculate the L_2 oriented distance between them by default. It will calculate the unoriented distance if orient is set to False (or both projected forms share a sign or either are achiral), and the oriented or unoriented Chebyshev distance if dtype is set to 0


# Coform
Calculation of, and with coforms and voforms from lattices. Requires standard packages, plus Pymatgen for the generation of Niggli Reduced cells and the general calculation of lattices from parameter. 

UTILITY FUNCTIONS:

The utility function *tolcheck* is used to check closeness of values to within a certain tolerance. 

The utility function *dist* takes any two lists of the same length dst(a,b) and a type 0, 1 or 2 
It calculates the L_inf, L_1 or L_2 distances between the lists as an array

The utility function *canon_coform* takes a root form (as a 3x2 array) and returns it in its canonical version (that is, sorted such that the leftmost two entries in the top column are the smallest, and ordered). 

The utility function *baryproj* projects any triple into barycentric co-ordinates in the quotient triangle. If all the input triples os ordered, the point will be in the qutient triangle. Similarly, *orthproj* make an orthogonal projection along three vectors at angles of 2pi/3. 

You shouldn't need to use *one_red_step* unless you want to run some sort of unit step - it just runs a single coform reduction. 

The utility function *lat_param* will take lattice parameters in the form (a,b,c, alpha, beta, gamma) and produce the vectors of the **niggli reduced** lattice in a format suitable to be taken as a Lattice_3D object

OBJECTS:

The *Lattice_3d* object should have the form of a list of d co-ordinates in R^3. The most important function is *get_lat_rootform*, which will return the root form of any lattice object in its canonical permutation (it will perform a Niggli reduction automatically first). If you just want the original coform  *get_lat_coform* will do the job. Hopefully every other function is self-explanatory. 

In both cases, the output is a *Coform* object (which if entered directly should be a 2 x 3 array). Its associated *Voform* can be calculated by the *get_voform* function. If you are entering a 'putative' non-reduced coform  the *red_cof* function will produce the actual positive coform in canonical form, while *cf_root* converts that into a root form and *cf_root_bary* will plot barycentric co-ordinates for the top row of the coform in the canonical. The *cf_permute* function outputs all permitted permutations of the coform values as a vector in R^6. 

The *Voform* object should be entered in the format [[v_0^2], [v_1^2, v_2^2. v_3^2], [v_01^2. v_02^2, v_03^2]]. The *get_coform* function calculates the associated coform.

A coform and voform together can be turned into a *CofVof* object. It can be used to reconstruct a lattice directly from coform and voform values: *get_angles*, *get_lat_lens* and *get_lattice_vectors* all do exactly what they say they do, the latter outputting a list that you can make into a Lattice_3d object if you like. The *full_reduction* function takes any **putative** coform-voform pair and, if not already reduced, returns the reduced pair. 


DISTANCE CALCULATION:
There are two coform distance calculator functions

*lat_cfdist* takes two Lattice_3d objects as input. dtype defaults to 2 and calculates the coform distance derived from the standard Euclidean distance. Note that it therefore needs to go through as many reduction steps as are required to get to positive conorm values. Set dtype = 1 for taxicab metric (L_1) and dtype = 0 for max metric (L_inf). 

*cfdist* does the same thing, but takes two Coform objects as input. 
