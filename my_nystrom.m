function [U, Sig, V] = my_nystrom(B, m, k, p)
if isnumeric(B)
    A =@(x) B * x;
else
    A =@(x) B(x);
end

l = p + k;
G = randn(m, l);
X = A(G);
C = G'*X;
[Q, R] = qr(X, 0);
M = R * (C\R'); M = (M + M')/2;
[u, s] = svd(M);
U = Q * u(:, 1 : k);
V = U;
Sig = s(1 : k, 1 : k);