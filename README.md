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

In both cases, the output is a *Coform* object (which if entered directly should be a 2 x 3 array). Its associated *Voform* can be calculated by the *get_voform* function. If you are entering a 'putative' non-reduced coform  the *red_cof* function will produce the actual positive coform in canonical form, while *cf_root* converts that into a root form and *cf_root_bary* will plot barycentric co-ordinates for the top row of the coform in the canonical 

The *Voform* object should be entered in the format [[v_0^2], [v_1^2, v_2^2. v_3^2], [v_01^2. v_02^2, v_03^2]]. The *get_coform* function calculates the associated coform.

A coform and voform together can be turned into a *CofVof* object. It can be used to reconstruct a lattice directly from coform and voform values: *get_angles*, *get_lat_lens* and *get_lattice_vectors* all do exactly what they say they do, the latter outputting a list that you can make into a Lattice_3d object if you like. The *full_reduction* function takes any **putative** coform-voform pair and, if not already reduced, returns the reduced pair. 


DISTANCE CALCULATION:
There are two coform distance calculator functions

*lat_cfdist* takes two Lattice_3d objects as input. dtype defaults to 2 and calculates the coform distance derived from the standard euclidean distance. Set dtype = 1 for taxicab metric (L_1) and dtype = 0 for max metric (L_inf). Set the parameter reflect = True if you wish to include reflections as part of the permutation group. 
