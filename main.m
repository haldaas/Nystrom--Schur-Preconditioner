clear all; clc;
format shortE;

GEVP = 1; % ideal two-level preconditioner
Nystrom = 1; % Nyström--Schur preconditioner
Jacobi = 1; % One-level preconditioner
hsl = 0; % Matlab interface for HSL_MI28 is required
iChol = 1; % Matlab incomplete Cholesky

k = 20; % Dimension of the deflation space (the rank of the correction term)
p = 0; % Oversampling parameter for Nyström's method

% params for the approximate solution of S_I system
mid_tol = 0.1;
mid_maxit = 500;

rootdir = './';
addpath(rootdir);

% solver
tol = 1e-6;    % convergence tolerance
MaxIt = 1500;  % maximum iteration count

% configuration
level = 6;     % 2^level is the number of subdomains (number of diagonal blocks)

% setup
nbSubdomain = 2^level;
N = nbSubdomain;
nbPartition = nbSubdomain + 1;

% read matrix:
%%%%%%%%%%%%%%%%%%%%%
%To change the matrix load the matrix required and store it in the variable AP
load s3rmt3m3.mat;
AP = Problem.A;
%%%%%%%%%%%%%%%%%%%%%

rng(1);
bP = randn(length(AP),1);

fprintf("Problem size: %d\n", length(AP));

%call 1-level K-way nested dissection
[permutation,partsize,partbegin,partend] = OneLevelKwayNestedDissection(AP, level);

%permute the matrix
A = AP(permutation, permutation);

b = bP(permutation);
[A, b, R, C] = equil_rar_Ab(A, b);
A = (A + A')/2;

n = size(A,1);

% set blocks
D = cell(nbPartition, 1);
RD= D;
for i = 1 : nbPartition
   D{i} = A(partbegin(i) : partend(i), partbegin(i) : partend(i));
   RD{i} = chol(D{i});
end

L = cell(nbSubdomain, 1);
for i = 1 : nbSubdomain
   L{i} = A(partbegin(nbPartition) : partend(nbPartition), partbegin(i) : partend(i));
end

% set Schur complement
Ag = D{nbPartition}; ng = length(Ag);
fprintf("Schur complement size: %d\n", ng);

RAg = chol(Ag);
Agm1 =@(x) solveAg_SPD(RAg, x);
S =@(x) applyS_SPD(Ag, RD, L, nbSubdomain, x);
AIm1 =@(x) solveAI_SPD(RD, N, partbegin, partend, x);
AI =@(x) applyAI_SPD(D, N, partbegin, partend, x);
AgI =@(x) applyL_SPD(L, partbegin, partend, N, x);
AIg =@(x) applyU_SPD(L, partbegin, partend, N, x);
BtDm1B =@(x) applyBtDm1B_SPD(RD, L, N, x);

bI = bP(1:n - ng);
bg = bP(n-ng + 1 : n);

% params for the approximate solution of S_I system
prec =@(x) AIm1(x);

op1 =@(x) AI(x) - AIg(Agm1(AgI(x)));
Op1 =@(x) AgI(myBFBCG(op1, AIg(x), mid_tol, mid_maxit, prec));
Op2 =@(x) Op1(x')';

my_legend = [];

if(GEVP)
    opts.disp = 1;
    opts.tol = 1e-4;
    opts.maxit = 600;
   [Z, Lam] = eigs(BtDm1B, ng, Ag, k+p, 'largestabs', opts);

    lam = diag(Lam);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    E = Lam * diag(1./(1 - lam));
    M1 = @(x) Agm1(x) + Z * (E*(Z'*x));
    Mm1 =@(x) applyMm1_SPD(L, AIm1, M1, partbegin, partend, N, x(1:n-ng), x(n-ng+1:n));
    
    [xG,flag,relres,iter_GEVP,resvecGEVP] = pcg(A, b, tol, MaxIt, Mm1);

    resvecGEVP = resvecGEVP/resvecGEVP(1);
    %%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%
    my_legend = [my_legend, "GEVP"];
end
%%%%%%%%%%%%%%%%%%%%%%
if(Nystrom)
    [Z, Sig, V] = my_nystrom(Op1, ng, k, p);
    Z = Agm1(Z);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    M1 = @(x) Agm1(x) + Z * (Sig*(Z'*x));
    Mm1 =@(x) applyMm1_SPD(L, AIm1, M1, partbegin, partend, N, x(1:n-ng), x(n-ng+1:n));
    
    [xN,flag,relres,iter_Nystrom,resvecNystrom] = pcg(A, b, tol, MaxIt, Mm1);
    
    resvecNystrom = resvecNystrom/resvecNystrom(1);
    my_legend = [my_legend, 'Nystr\"om--Schur'];
end
%%%%%%%%%%%%%%%%%%%%%%
if(Jacobi)
    M1 =@(x) Agm1(x);
    Mm1 =@(x) applyMm1_SPD(L, AIm1, M1, partbegin, partend, N, x(1:n-ng,:), x(n-ng+1:n, :));
    [xBJ,flag,relres,iter_BJacboi,resvecBJacobi] = pcg(A, b, tol, MaxIt, Mm1);
    resvecBJacobi = resvecBJacobi/resvecBJacobi(1);
    my_legend = [my_legend, "BJacobi"];
end
%%%%%%%%%%%%%%%%%%%%%%
if(hsl == 1)
    pciL = hsl_mi28_precond(AP);
    [xho,flag,relres,iter_hslo,resvechslo] = pcg(AP, bP, tol, MaxIt, @(x) pciL.apply(x));
    resvechslo = resvechslo/resvechslo(1);
    my_legend = [my_legend, "\verb|HSL_MI28|"];
end
%%%%%%%%%%%%%%%%%%%%%%
if(iChol == 1)
    iL = ichol(AP,struct('diagcomp',.1));
    [xc,flag,relres,iter_ichol,resvecichol] = pcg(AP, bP, tol, MaxIt, iL, iL');
    resvecichol = resvecichol/resvecichol(1);
    my_legend = [my_legend, "\verb|ichol|"];
end

% plotting
figure
iters = [];
residuals = [];
itersName = {};
normb= norm(b);
if(GEVP)
   semilogy(resvecGEVP, '-b');
    hold on; grid on;
    iters = [iters, iter_GEVP];
    itersName = [itersName, 'GEVP'];
    residuals = [residuals, norm(b-A*xG)/normb];
end
if(Nystrom)
    semilogy(resvecNystrom, '-k');
    hold on; grid on;
    iters = [iters, iter_Nystrom];
    itersName = [itersName, 'Nystrom'];
    residuals = [residuals, norm(b-A*xN)/normb];
end
if(Jacobi)
    semilogy(resvecBJacobi, '--kd');
    iters = [iters, iter_BJacboi];
    itersName = [itersName, 'BJacobi'];
    residuals = [residuals, norm(b-A*xBJ)/normb];
end
if(hsl)
    semilogy(resvechslo, ':k');
    iters = [iters, iter_hslo];
    itersName = [itersName, 'hsl'];
    residuals = [residuals, norm(bP-AP*xho)/norm(bP)];
end
if(iChol)
    semilogy(resvecichol, '--k');
    iters = [iters, iter_ichol];
    itersName = [itersName, 'iChol'];
    residuals = [residuals, norm(bP-AP*xc)/norm(bP)];
end

legend(my_legend,'Interpreter','latex')
for i = 1 : length(iters)
    fprintf("%s\t\t", itersName{i});
end
fprintf("\n");
for i = 1 : length(iters)
    fprintf("%d \t\t", iters(i));
end
fprintf("\n");
for i = 1 : length(iters)
    fprintf("%.1e \t", residuals(i));
end
fprintf("\n");
