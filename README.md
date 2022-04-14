# Rootform_2d
Calculation of and with root forms in 2 Dimensions

**Utility Functions**

*distgen(q, v_1, v_2)* is a generic Minkowkski distance calculator in R^2. It takes two inputs - a value of q >= 0 and two points in R^2 given as lists. If q=0 it will return the Chebyshev (L_inf) distance between points, otherwise it will return the L_q distance. 

*mincyc(l)* takes a list l of any length and returns the cyclic permutation that puts the smallest list member first

*roundlist(l,r)* takes a list of floats l and an integer value r. It returns l with r-rounded float values. 

*minf(a,b,c,d)* calculates the quantity max(|a-b|, |c-d|, |a+b-c+d|/2) required for the root form distance

*makelat(a,b,t)* will return the vectors of a Niggli Reduced 2D lattice given two lengths a, b and an angle t

*sb_sign(veclist)* will take a list of three vectors and, assuming that they are an obtuse superbase, calculate the sign of that superbase

*index_sorted(list)* will take any list of quantities, and return a related list whose value at each index is the position the quantity at that index in the original list would be at in a sorted version of the list. 

**Lat2D**

The Lat2D class should be entered as two lists [[x_1, y_1], [x_2, y_2]]. The following methods are available:

*make_obsb()* will iterate the reduction step detailed in the paper until an obtuse superbase is reached

*sb_sign()* will calculate the sign of the superbase (positive if superbase vectors are ordered anticlockwise, negative if superbase vectors are ordered clockwise)

*make_cf* returns the **unsorted** coform for the lattice

*lattice_sign(tol = 10**-6)* returns the sign of the lattice as calculated from the coform - it reverses the sign of the superbase if a swap of two values is required to order the coform values list, or keeps that sign if ordering the list requires only permutation. It will the sign to 0 if any value in the coform is less than the input tolerance (10**-6 default). 

*make_rf(tol= 10**-6* returns the signed root form object - it will return 0 if the coform is achiral (to within a set tolerance of $10^-6$), otherwise it will return the sign passed to it by the *lattice_sign* function. 

**RF_2signed Class**

RF2_signed objects consist of a pair of inputs - an ordered list of three positive numbers and a sign (0, +1, -1). If an **achiral** root forms is entered, the sign will default to 0 regardless of the input. If the user inputs an unordered list there will be a warning message:

*rightsign()* will return the unordered root form - that is, it will swap the last two values of the root form if the sign is negative

*rf_chir(dtype = 0, pgroup = 2)* will, return the signed L_inf or L_2 based point group chirality of the root form. Id pgroup is set to 4 or 6, the D4 or D6 chirality will instead be calculated. If dtype is set to 2, the L_2 distance will be calculated rather than L_inf

*projform()* will return the co-ordinates of the form in the quotient triangle, along with the sign of the root form, as a PF2 object. 

*coform2d()* and *voform2d()* will return the coform  and voform for the obtuse superbase represented by the root form. Vectors are in anticlockwise order as stated in the paper. 

*make2lat()* will reconstitute the actual 2 dimensional lattice based on the input root form. Note that the output lattice will **not** be Niggli reduced - it will be the two vectors v_1, v_2 calculated as part of the obtuse superbase. 

**PF2 Class**

Projected Form (PF2) objects should be entered  as a list [x, y] followed by a value +/-1 or 0 indicating the sign of the lattice from which the projected form is derived, and a tolerance set by defalut at 10**-6 - the sign will set to zero if either x or y value is within the tolerance level. The following methods are available:

*qs_plot()* will return the co-ordinates of the projected form in the quotient square. 

*pf_grpchir(dtype = 0, pgroup = 2)* returns the signed L_inf or L_2 chirality measure for the projected form depending on whether dtype is set to 0 or 2. Setting pgroup to 4 or 6 will calculate the D4 or D6 Chirality. Any other setting of dtype or pgroup will return an error message. 

*pf_chir(dtype = 0)* will return the signed D2 chirality. Defaults to use of Chebyshev distance - if dtype = 2 the Euclidean distane will be used instead. 

*root_from_pf2 (scale = 1)* will return a signed root form derived from an input projected form, which can be scaled by setting the scale factor

*lattice_from_pf2 (sc=1)* will return the 2D lattice reconstituted from an input projected form - setting the scale factor scales the underlying root form

*sphere_proj()* returns the latitude and longitude of the spherical projection. 

**Distance Calculations**

*pf2dist(pf_1, pf_2, orient = True, dtype = 2)* will take two projected forms and calculate the L_2 oriented distance between them by default. It will calculate the unoriented distance if orient is set to False (or both projected forms share a sign or either are achiral), and the oriented or unoriented Chebyshev distance if dtype is set to 0

*rf2dist(rf_1, rf_2, orient = True, dtype = 2)* will take two root forms and calculate the L_2 oriented distance between them by default. It will calculate the unoriented distance if orient is set to False (or both projected forms share a sign or either are achiral), and the oriented or unoriented Chebyshev distance if dtype is set to 0

