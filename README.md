# Lattice_Invariance
Calculation of Conorm and Vonorm Invariants

The utility function *dist* takes any two lists of the same length dst(a,b) and a type 0, 1 or 2 
It calculates the L_inf, L_1 or L_2 distances between the lists as an array

The utility function *rebalance* takes a list, an index and a value. It subtracts the given value from the list element at the given index,
and adds that value to the elements at all other indices

The utility function *swapper* takes a list and two elements and swaps the elements in the list

The utility functions *angcalc_2d* and *angcalc_3d* [TBD] take lattice vector lengths and angles and output the lattice as a triple of cartesian vectors, with the first vector aligned to the x axis and (for 3d lattices) the first and second axes spanning the xy plane. Vectors default to unit length. 

*Lattice_2d* and *Lattice_3d* are lattice objects. Theytake a list of d co-ordintes in R^d, eithe directly or as the output of the angcalc utility functions.
The method msb/msb3d makes a superbase out of a 2D or 3D lattice simply by adding a v_0, which is the negative sum of all the input vectors

*Obsuperbase_3d* and *Obsuperbase_2d* take what they assume to be an obtuse superbase and calculate Selling parameters

Their respective reduction methods *red2d* and *reduce_3d* run the algorithm decsribed in Conway and Sloane until all Selling parameters are positive. 

The methods *allperms_2d* and *allperms_3d* give a list of all permitted coform permutations for an input superbase (which does not necessarily have to be reduced). If the boolean input variable *iso* is true, coform permutations induced by lattice reflections are suppressed. 

The *isodist_2d* and *isodist_3d* functions take two lattice objects of the relevant dimension and a dtype variable of values 0,1,2 and a boolean variable *iso* as described above. It calculates the L_(dtype) distance between the two lattices (with reflections suppressed if iso = True). If dtype = 0, the L_infinity distance is calculated.

