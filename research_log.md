#### May 06,2023

1. Write efficient code
2. Figure out how to go from HOI RE to HOI GLV in a systematic way
3. run 1000 simulations with d = 2,3,4,5,6 and n = 2,4,8,16 on the cluster

#### May 08,2023

1. Find roots in mathematica too, in order to check that our results are correct
2. Find a way to prove random zoo with three-way interactions. Understand Carlos's proof
3. Numerically get all real equilibria, multiply by all possible D's, and see what is the probability of feasibility

#### May 09,2023

test if my results change with symmetric tensor slices.
Test for dimension 2, different variances
build output name outomatically

#### May 10,2023

Compare results with d = 2 to see the increase, and base research in known fact(s (1/2^n)
Check variance also converges to the theoretical formula
Put errorbars in simulations
Do the same for Uniform distribution
Write down a way go from RE to GLV for arbitrarily high order interactions

#### May 11,2023

Show that GLV with HOIs is a random polynomial
Compute probability of feasibility analytically (use Jensen's inequality, maybe ping Alex)
Show in simulations the probability of feasibility

#### May 12,2023

heat map with probability of feasibility as a function of d and n. 

#### May 15,2023

code jacobian in julia to check for stability after each simulation is finished. 
To do this, first check that polynomial way of writing this works, then work with the polynomials, by taking the derivative in Julia
the name of the file with stability data is "dim_6_div_6_s_1000_normal"
find a way of, given r, A, B, ..., convert to a polynomial in julia

#### Jul 20,2023

figure out why certifysolutions is not working properly (case of i = 3)

#### Aug 16,2023

    Plan for new set of simulations: 
    1. Calculate jacobian of system
    2. Evaluate at all feasible equilibria
    3. Get dominant eigenvalue
    4. If its negative, classify as stable
    4. Save equilibria
    5. If its possitive, run dynamics and record the equilibrium that system reaches. 
    6. Trim the system, and calculate all feasible equilibria there, and thier local stability
    7. Save equilibria
    8. Save all parameters used as files capable to be loaded again

    There are then three type of files: 
    1. Summary statistics: n, d, nsol, npos
    2. equilibria values: n, d, eq1, domeig1
    		          n, d, eq2, domeig2
			  ...
    3. parameter values: save big tensor as a julia file, with seed name, to be able to load again
