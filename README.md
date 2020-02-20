IntervalAssignment 
GitHub repository
https://github.com/samitch/IntervalAssignment.git 

IIA, Incremental Interval Assignment by Integer Linear Algebra, is a solver for optimizing integer matrix problems Ax=b.
IIA was developed for the application of deciding the number of mesh edges on model curves (intervals) for quad and hex meshing. Meshing schemes impose constraints, ranging from the mild requirement that any quad mesh must have an even number of edges on its boundary, to structured mapped patches where opposite sides of a rectangle must have exactly equal numbers of edges.


== Problem & Solution

The Interval Assignment (IA) problem that IIA solves is

min lex f(x), where f(x) = (x>g ? x/g : g/x)
s.t. Ax = b
     x in [lo,hi]
     A, x, b integer

Inputs are matrix A, and vectors b, g, lo, hi, and for each row whether "=" in Ax=b is equality or inequality or "sum-is-even".
Outputs are vector x.

Here f(x_i) is the ratio of x_i to g_i. Here g is the goal, the user-defined ideal number of intervals for each variable, and may be positive floating point or "don't care," in which case that variable does not contribute to f(x). 

Note all variables may be set to "don't care," in which case IIA just solves the constraints. It may be possible to use IIA to solve other optimization problems of IA's form that do not arise from the context of interval assignment, but this is untested. Features of the IA context are A is sparse and A and b entries are small integers: A non-zeros are mostly +/-1 and a few are +/-2, and b is mostly 0. Also lo is typically 1 and hi is typically unbounded, except for "don't care" and slack variables. It may be possible to edit the source code to use a different function for f, but this is untested. The algorithm implementation is fundamentally based on the objective being min lex, the minimum lexicographic vector of f.

IIA is a discrete algorithm based on sparse integer linear algebra, variants of Gaussian elimination. The only floating point computations involve g, which is used in defining downhill directions and the lexicographic order.

IIA process outline
0. convert <= and "sum-is-even" constraints to equality = using slack variables. (robust)
1. solve Ax = b using Hermite Normal Form (HNF) to find x1. (robust)
2. solve Ax = 0 using Reduced Row Echelon Form (RREF) to find vectors N spanning the nullspace of A. (robust)
3. solve x in [lo,hi] by adding linear combinations of N to x1, in downhill directions, to find x3. (can fail)
4. min lex f(x) by adding linear combinations of N to x3, in downhill directions constrained by [lo,hi], to output x4. (can fail)
Steps 1-4 are performed over (semi) independent subproblems and recombined for the global solution.

The main source of non-robustness is finding linear combinations of N that point in downhill directions. The heuristic is based on Gaussian elimination to avoid uphill directions. We save uphill directions found in previous searches, so the search is not exhaustive. This keeps running time low (sub-exponential) but can cause the algorithm to terminate with a sub-optimal solution, and possibly one that does not satisfy x in [lo,hi].


== Why Use It?

The benefits of IIA are speed and always producing an integer solution. The drawback is the potential for the bounds to be unsatisfied or the solution to be sub-optimal. Example serial runtime is 0.2 seconds for matrix A sized 2000 x 2000 with about 4000 non-zeros. This is phenomenally fast for an Ax=b integer problem. In contrast, a predecessor code (BBIA) based on linear programming and branch and bound solves the same example problem in 20 minutes, 6000x slower.

IIA is C++ and has *no* dependencies. 
IA.h is the interface.
The driver code test.cpp gives examples of setting up and solving the problem.


== Acknowledgements

The predecessors to this code include BBIA, a linear-programming plus branch-and-bound solver, implemented in Cubit in 1995--1997 by Scott A. Mitchell and maintained ever since (2020+) by the Cubit team [BBIA]. Another predecessor code by Scott A. Mitchell is NLIA, a non-linear programming based approach that attempted to reduce (but not eliminate) the need for branch and bound [NLIA]. A direct progenitor is MSIIA, a version developed for the simpler context of refining or coarsening an existing hex mesh in a semi-structured way [MSIIA, MSIIA LDRD]. This code (IIA) was developed in Cubit, then extracted to become this library-like standalone.

"I" is the IIA author, Scott A. Mitchell. I thank Jason Shepherd, Robert Kerr, Michael Plooster, and Clinton Stimpson for their work on BBIA, and the related work of defining what constraints and goals should be sent to the solver. I thank Clinton Stimpson for his work on the infrastructure related to the Cubit version of IIA. I thank Timothy Tautges for prioritizing the interval assignment problem and supporting the work of both BBIA and NLIA. I thank Matthew Staten for supporting and inspiring MSIIA, which was the genesis for me thinking IIA might be a viable approach for the general problem. I thank Roshan Quadros, Michael Parks, and Michael Skroch for suporting the development of IIA. I thank Roshan Quadros, Trevor Hensley, and Salome Thorson for faciliating discussions with Cubit users  Jacquelyn Rae Moore, Neal Grieb, and others, to help determine desirable constraints, bounds and options. I thank David White for his work on autoscheme, the problem of deciding which meshing scheme (algorithm) to use for individual surfaces and volumes, and for using interval assignment to determine if a scheme is feasible. I thank Paul Stallings and Byron Hanks for discussion on how to define a clean interface to IIA. To the best of my knowledge there are no other authors of IIA besides Scott A. Mitchell, with the exception of some generic infrastructure from Cubit developed by the people mentioned in this acknowledgement.


== Bibliography

[BBIA]
  "High Fidelity Interval Assignment," Scott A. Mitchell, Proc. 6th International Meshing Roundtable `97, 33-44 (1997), and
International Journal of Computational Geometry and Applications, Vol. 10, No. 4 (2000) 399-415.

[NLIA] 
  "Simple and Fast Interval Assignment Using Nonlinear and Piecewise Linear Objectives," Scott A. Mitchell, Proc. 28th International Meshing Roundtable, year 2013, vol 22, pages 203-221, eds. Josep Sarrate and Matthew Staten, Springer, ISBN 978-3-319-02335-9, doi 10.1007/978-3-319-02335-9, url http://www.cs.sandia.gov/~samitch/bibliography_2007.html.

[MSIIA]
  "Incremental Interval Assignment for Mesh Scaling Assembly Models," research abstract, Scott A. Mitchell, International Meshing Roundtable, 2019. https://doi.org/10.5281/zenodo.3653101

[MSIIA LDRD]
  "Incremental Interval Assignment (IIA) for Scalable Mesh Preparation," Scott A. Mitchell, Sandia National Laboratories technical report SAND2019-11163R, LDRD 19-1050, 2019.

[IIA]
  "Interval Assignment by Integer Linear Algebra," Scott A. Mitchell, in preparation.


== IIA Design Background 

The following online lecture notes and discussions were helpful in designing IIA 

HNF Hermite Normal Form and why it is what we want to solve Ax=b
  http://sites.math.rutgers.edu/~sk1233/courses/ANT-F14/lec3.pdf
  http://sites.math.rutgers.edu/~sk1233/courses/ANT-F14/lec4.pdf

RREF Reduced Row Echelon Form
  https://people.sc.fsu.edu/~jburkardt/c_src/row_echelon_integer/row_echelon_integer.html

How RREF helps us compute the nullspace of a maxtrix
  https://math.stackexchange.com/questions/88301/finding-the-basis-of-a-null-space

How the nullspace helps us solve integer optimization problems
  https://cs.stackexchange.com/questions/82086/linear-programming-restricted-to-rational-coefficients

  Slide 27, totally unimodular ILP solve (IA's matrix A is not totally unimodular, but is sometimes close)
  https://www.isical.ac.in/~arijit/courses/autumn2016/ILP-Lecture-1.pdf
