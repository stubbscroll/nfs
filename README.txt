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
- Algebraic square root using the algorithm by Couveignes which requires that
  f(x) has odd degree

Read more about the number field sieve elsewhere (for example
http://en.wikipedia.org/wiki/General_number_field_sieve ). Be warned, the
algorithm is extremely complicated and math-heavy.

It is made to be dependent on as few external libraries as possible.
Currently, the only requirement is the GMP (GNU Multiprecision) library.
The programming language is C, and the only C99 feature used is long long.
It means that gcc will happily compile the code as C89.

The program is fully working, but with the following limitations:
- f(x) must have degree 3. The polynomial factorization algorithm only
  searches for linear factors, and a general routine is needed before
  degrees of 4 or larger can work.
- f(x) must have odd degree. This requirement is because of the algebraic
  square root algorithm.

Bottlenecks:
- Linear algebra step (Gauss-Jordan for huge matrices is slow).
- Algebraic square root. Some shortcuts were taken when implementing the
  algorithm, which made the algorithm slower.
- Hence, the program will not work well for n larger than 50-60 digits.

Areas of improvement:
- Profile the program and optimize the hot spots
- Factorize the constant term of f(x) and use this factorization to generate
  all divisors. In this way make a faster check for linear factors that should
  work for larger numbers (this improvement is made superfluous by proper
  polynomial factorization)
- Use large primes in the sieve to obtain higher yield (large prime variation)
- Use proper polynomial factorization of f(x), either with the extremely
  complicated Algorithm 3.5.7 in [Coh93], the even more complicated
  algorithm of [Len82] (which is said to be much slower than Algorithm 3.5.7)
  or use the naive, slightly dubious but easiest-to-implement algorithm
  described in section 3.5.5 in [Coh93].
- Use a faster linear algebra algorithm than Gauss-Jordan (for instance
  Block Lanczos or Block Wiedemann)
- Use lattice sieve instead of line sieve (expected speedup: 2-3x)
- Use a different algebraic square root algorithm that handles even degrees
- Better selection of f(x)
- Postprocessing of relations before linear algebra step

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
