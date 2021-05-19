# Lattice_Invariance
Calculation of Conorm and Vonorm Invariants

The utility function takes any two lists of the same length dst(a,b) and a type 0, 1 or 2 
It calculates the L_inf, L_1 or L_2 distances between the lists as an array

The utility function minperdist again takes two lists of the same length and type 0,1 or 2
It returns the minimum L_inf, L_1 or L_2 between all permutations of the first list and the second

Lattice_2d and Lattice_3d are lattice objects. Theytake a list of three Cartesian Co-ordinates. 
The method msb/msb3d makes a superbase out of it simply by adding a v_0, which is the negative sum of all the input vectors

Obsuperbase_3d and Obsuperbase_2d take what they assume to be an obtuse superbase and calculate Selling parameters
