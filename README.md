# Nystrom--Schur-Preconditioner

## __IMPORTANT:__
This code is for reproducing some numerical experiments in [1] and __SHOULD NOT__ be used otherwise.

## Introduction
This is a proof of concept Matlab code for the two-level Nyström--Schur preconditioner (NS-prec) proposed in [1]

The code requires a Matlab interface for metis 4. Download [metis 4](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/OLD/metis-4.0.3.tar.gz) and [metismex.c](https://www.cs.ubc.ca/~rbridson/download/metismex.c) then compile it on your platform.

To run comparisons against [HSL_MI28](https://www.hsl.rl.ac.uk/catalogue/hsl_mi28.html), the [Matlab interface for HSL_MI28](https://www.hsl.rl.ac.uk/catalogue/hsl_mi28.html) is required (a free academic licence can be obtained; see [licence](https://www.hsl.rl.ac.uk/download/HSL_MI28/2.2.2/)).

## Running test example
The file main.m is a script that compares the two-level NS-prec against the following preconditioners:
- the ideal two-level preconditioner
- the one-level preconditioner
- [HSL_MI28](https://www.hsl.rl.ac.uk/catalogue/hsl_mi28.html) (incomplete Cholesky factorization) -- requires Matlab interface; see [HSL_MI28](https://www.hsl.rl.ac.uk/catalogue/hsl_mi28.html)
- Matlab incomplete Cholesky factorization

One or multiple preconditioners can be included in the comparison by setting their value to 1 in the beginning of the file main.m

The parameters of NS-prec are the following:
- level: 2^level is the number of subdomains
- k: Dimension of the deflation space (20 is usually enough)
- p: Oversampling parameters for the Nyström method (0 is usually enough)

The repository includes the test matrix s3rmt3m3. To run the comparison on another matrix, follow the instructions given in the read matrix section inside the main.m file.

    [1] Two-level Nyström--Schur preconditioner for sparse symmetric positive definite matrices. 2021.
        Hussam Al Daas and Tyrone Rees and Jennifer Scott
    
 To cite this work:
 
    @article{doi:10.1137/21M139548X,
    author = {Al Daas, Hussam and Rees, Tyrone and Scott, Jennifer},
    title = {Two-Level Nyström--Schur Preconditioner for Sparse Symmetric Positive Definite Matrices},
    journal = {SIAM Journal on Scientific Computing},
    volume = {43},
    number = {6},
    pages = {A3837--A3861},
    year = {2021},
    doi = {10.1137/21M139548X},
    URL = {https://doi.org/10.1137/21M139548X},
    eprint = {https://doi.org/10.1137/21M139548X}
    }
