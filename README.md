# Nystrom--Schur-Preconditioner

## __IMPORTANT:__
This code is for reproducing some numerical experiments in [1] and should not be used otherwise.

This is a proof of concept Matlab code for the two-level Nyström--Schur preconditioner (NS-prec) proposed in [1]

The code works on Linux and mac platforms. To run on windows, a Matlab interface for metis 4 on windows is required (this is not tested).
To run comparisons against [HSL_MI28](https://www.hsl.rl.ac.uk/catalogue/hsl_mi28.html), the [Matlab interface for HSL_MI28](https://www.hsl.rl.ac.uk/catalogue/hsl_mi28.html) is required (a free license for personal use).

## Running test example
The file main.m is a script that compares the two-level NS-prec against the following preconditioners:
- the ideal two-level preconditioner
- the one-level preconditioner
- [HSL_MI28](https://www.hsl.rl.ac.uk/catalogue/hsl_mi28.html) (incomplete Cholesky factorization) -- requires Matlab interface; see [HSL_MI28](https://www.hsl.rl.ac.uk/catalogue/hsl_mi28.html)
- Matlab incomplete Cholesky factorization

The parameters of NS-prec are the following:
- level: 2^level is the number of subdomains
- k: Dimension of the deflation space (20 is usually enough)
- p: Oversampling parameters for the Nyström method (0 is usually enough)

The repository includes the test matrix s3rmt3m3. To run the comparison on another matrix, follow the instructions given in the read matrix section inside the main.m file.

    [1] Two-level Nyström--Schur preconditioner for sparse symmetric positive definite matrices. 2021.
        Hussam Al Daas and Tyrone Rees and Jennifer Scott
    
 To cite this work:

    @misc{daas2021twolevel,
      title={Two-level Nystr\"om--Schur preconditioner for sparse symmetric positive definite matrices}, 
      author={Hussam Al Daas and Tyrone Rees and Jennifer Scott},
      year={2021},
      eprint={2101.12164},
      archivePrefix={arXiv},
      primaryClass={math.NA}}
