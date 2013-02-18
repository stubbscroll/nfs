Number field sieve v0.0 by Ruben Spaans

This is a simple implementation of the basic variant of the number field sieve
algorithm for factorizing integers. "Basic" implies that the easiest algorithm
is usually chosen at each step. The program currently consists of the following
sub-algorithms:

- Find minimal polynomial f(x) using base-m expansion of n
- Factorization of f(x) that searches for linear factors only; runs only if the
  constant term is less than 2*10^9. Limits f(x) to degree 3
- Find rational factor base using Sieve of Eratosthenes
- Find pairs (p,r) for the algebraic factor base by using an efficient root
  finding algorithm for polynomials over Z_p
- Line sieve using approximations of base 2 logarithms
- Solve matrix step using Gauss-Jordan
- Algebraic square root is currently not implemented

It is made to be dependent of as few external libraries as possible.
Currently, the only requirement is the GMP (GNU Multiprecision) library.

Things needed for a working program:
- Polynomial factorization of f(x) that works in all cases (currently only
  works if constant term is less than 2*10^9 and finds linear factors only)
- Algebraic square root

Limitations:
- Currently limited to degree 3 polynomials, because the factorization of
  f(x) only looks for linear factors.
- Gauss-Jordan is a memory bottleneck, limiting the size of the integers to
  factor to around 60 digits on a machine with 16 GB RAM

Areas of improvement:
- Profile the program and optimize the hot spots
- Use large primes in the sieve to obtain higher yield
- Use proper polynomial factorization of f(x), either with the extremely
  complicated Algorithm 3.5.7 in [Coh93], the even more complicated
  algorithm of [Len82] (which is said to be much slower than Algorithm 3.5.7)
  or use the naive, slightly dubious but easiest-to-implement algorithm
  described in section 3.5.5 in [Coh93].
- Use a faster linear algebra algorithm than Gauss-Jordan (for instance
  Block Lanczos or Block Wiedemann)
- Use lattice sieve instead of line sieve
- Use a different algebraic square root algorithm that handles even degree
- Better selection of f(x)

Usage:

nfs < inputfile

where the inputfile is a text file containing the following information on
separate lines:

- Integer to factorize
- Upper bound of algebraic factor base (or 0 to let the program decide)
- Upper bound of rational factor base (or 0 to let the program decide)
- Number of quadratic characters (or 0 to let the program decide)
- Degree of polynomial (odd number)
- m value for base-m expansion in polynomial creation (or 0 to let
  the program decide). m much result in a monic polynomial of the degree as
  specified above.
- Sieve width a (sieves from -a to a)
- Threshold for accepting numbers in the sieve (base 2 logarithm)
- Number of smallest primes to skip on each side of the sieve
- Number of extra relations wanted for linear algebra

The integer to be factorized can be specified literally, or by c[d] (a
random composite number with d digits and no very small factors) or by r[d]
(a random composite number that is the product of two similarly sized primes).

References:

[Coh93] A course in computational algebraic number theory (Henri Cohen, 1993)
[Len82] Factoring polynomials with rational coefficients (Lenstra, Lenstra,
        Lovacz, 1982)
