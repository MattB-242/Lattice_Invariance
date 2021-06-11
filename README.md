# Lattice_Invariance
Calculation of Conorm and Vonorm Invariants in 2 and 3 dimensions

UTILITY FUNCTIONS:

The utility function *dist* takes any two lists of the same length dst(a,b) and a type 0, 1 or 2 
It calculates the L_inf, L_1 or L_2 distances between the lists as an array

The utility function *rebalance* takes a list, an index and a value. It subtracts the given value from the list element at the given index,
and adds that value to the elements at all other indices

The utility function *swapper* takes a list and two elements and swaps the elements in the list

For the 2D code, the function *angcalc_2d*  take lattice vector lengths and angles and output the lattice as a triple of cartesian vectors, with the first vector aligned to the x axis and (for 3d lattices) the first and second axes spanning the xy plane. Vectors default to unit length. 

For the 3D code, the function *cartcalc_3d* encodes a library of Bravais lattices. It will output primitive unit cell vectors for a given lattice (cubic, rhombic, triclinic, hexagonal, orthorhombic, isoclinic) and an appropriate type (Primitive, Face-Centered, Base-Centered, Body_centered') with correctly input vector lengths and angles. 

*Lattice_2d* and *Lattice_3d* are lattice objects. Theytake a list of d co-ordintes in R^d, eithe directly or as the output of the angcalc utility functions. They both have superbase methods.

*Obsuperbase_3d* and *Obsuperbase_2d* take what they assume to be an obtuse superbase and calculate Selling parameters
Their respective reduction methods *red2d* and *reduce_3d* run the algorithm decsribed in Conway and Sloane until all Selling parameters are positive. 

The methods *allperms_2d* and *allperms_3d* give a list of all permitted coform permutations for an input superbase (which does not necessarily have to be reduced). If the boolean input variable *iso* is true, coform permutations induced by lattice reflections are suppressed. 

The *isodist_2d* and *isodist_3d* functions take two lattice objects of the relevant dimension and a dtype variable of values 0,1,2 and a boolean variable *iso* as described above. It calculates the L_(dtype) distance between the two lattices (with reflections suppressed if iso = True). If dtype = 0, the L_infinity distance is calculated.

