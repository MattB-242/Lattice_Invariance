# Rootform_2d
Calculation of and with root forms in 2 Dimensions

**Utility Functions**

*distgen* is a generic Minkowkski distance calculator in R^2. It takes two inputs - a value of q >= 0 and two points in R^2 given as lists. If q=0 it will return the Chebyshev (L_inf) distance between points, otherwise it will return the L_q distance. 

*mincyc* takes lists of any length and returns the cyclic permutation that puts the smallest list member first

*roundlist* takes a list of floats l and an integer value r. It returns l with r-rounded float values. 

*redstep* performs a single reduction step for an input coform (three member list). 

**Working with Root Forms**

Root Form objects consist of a pair of inputs - an ordered list of three positive numbers and a sign (+1, or -1). An unoriented root form should be entered with a positive sign. **Achiral** root forms should also be entered with a positive sign, although if this is not adhered to the sign will be corrected. The following methods can be applied to a root form:

*rightsign* will return the unordered root form - that is, it will swap the last two values of the root form if the sign is negative

*rf_chirality* will, return the signed L_inf or L_2 based chirality of the root form, depending on whether dtype is set to 0 or 2

*projform* will calculate the co-ordinates of the projective form in the glued quotient triangles of the paper

*coform2d* and *voform2d* will return the **oriented** coform and voform for the obtuse superbase represented by the root form. Vectors are in anticlockwise order as stated in the paper. 

*make2lat* will reconstitute the actual 2 dimensional lattice based on the input root form. 

**Working with Projected Forms**

Projected Form objects should be entered simply as a list [x, y]. The following methods are available:

*pf_chirality* returns the signed L_inf or L_2 chirality measure for the projected form depending on whether dtype is set to 0 or 2. Any other setting of dtype will return an error message. 

*root_from_pf2* will return the signed root form derived from an input projected form

*lattice_from_pf2* will return the 2D latice reconstituted from an input projected form 

[The following Root Form methods are planned for future updates]

*find_nearest_achiral_pf* will take a value q return the point that is L_q - closest to the boundary of the quotient triangle. If q=0 the Chebyshev distance will be returned.  

**Working with Lattices**

The Lat2D class should be entered as two lists [[x_1, y_1], [x_2, y_2]]. The following methods are available:

*make_CF* returns an anticlockwise oriented coform for the lattice. 

*make_RF* returns the signed root form. 


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
