function [x, varargout] = myBFBCG(Op, b, tol, maxit, prec, varargin)
% determin the type of operator and preconditioner
if(isa(Op,'float'))
    A =@(v) Op * v;
else
    A =@(v) Op(v);
end
if(isa(prec,'float'))
    M =@(v) prec \ v;
else
    M =@(v) prec(v);
end
if nargin == 6
    x0 = varargin{6};
    r = b - A(x0);
    counterOp = 1;
else
    x0 = zeros(size(b));
    r = b;
    counterOp = 0; 
end

% norm of rhs
normb = norm(b);
normr = norm(r);

% argout
flag = 4;
relres = normr/normb;
iter = 1;
resvec = zeros(maxit,1);
resvec(1) = normr;
blockSize = resvec;
blockSize(1) = size(b, 2);

% initial guess
x = x0;

% best approximate solution
bestx = x;
bestIter = 0;
bestres = normb;

% relative tolerance
reltol = normb * tol;
btol = normb * 1e-4;

k = 0;
Mr = M(r);
counterM = 1;
rtMr = r' * Mr;
p = Mr;
[p, s, ~] = svd(p, 0);
for j = 1 : size(s)
    if s(j, j) < btol
        fprintf("breakdown occured\n");
        break;
    end
end

blockSize(1) = j;
p = p(:, 1 : j);

while normr > reltol && k < maxit
    Ap = A(p);
    counterOp = counterOp + 1;
    temp = p' * Ap;
    alpha = temp\( p' * r);
    x = x + p * alpha;
    r = r - Ap * alpha;
    [~,sig,~] = svd(r,'econ');
%     diag(sig)'
    normr = norm(r);
    resvec(k + 2) = normr;
    if normr < resvec(bestIter + 1)
        bestx = x;
        bestIter = k + 1;
        bestres = normr;
    end
    Mr = M(r);
    counterM = counterM + 1;
    beta = - temp \ (Ap' * Mr);
    p = Mr + p * beta;
    [p, s, ~] = svd(p, 0);
    for j = 1 : size(s)
        if s(j, j) < btol
            fprintf("breakdown occured\n");
            break;
        end
    end
    blockSize(k + 2) = j;
    p = p(:, 1 : j);
    k = k + 1;
end

iter = k;
relres = resvec(bestIter + 1)/normb;
if(bestIter ~= k)
    fprintf("Best approximation given at iteration %d with relres %e\n", bestIter, bestres/normb);
    x = bestx;
end
if bestres < reltol
    flag = 0;
else
    flag = 1;
end
fprintf("iteration count %d with relres %e\n", iter, bestres/normb);

resvec = resvec(1 : iter + 1);
blockSize = blockSize(1 : iter + 1);
% blockSize'
if nargout > 1
    varargout{1} = flag;
end
if nargout > 2
    varargout{2} = relres;
end
if nargout > 3
    varargout{3} = iter;
end
if nargout > 4
    varargout{4} = resvec(1 : iter + 1);
end
if nargout > 5
    varargout{5} = blockSize(1 : iter + 1);
end
if nargout > 6
    varargout{6} = counterOp;
end
if nargout > 7
    varargout{7} = counterM;
end

% semilogy(resvec);
% plot(blockSize, 'o')
% hold on; grid on;