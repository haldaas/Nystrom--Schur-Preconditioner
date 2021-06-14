function [eqA,eqb, R,C] = equil_rar_Ab (A, b)

thresh = 0.1;
[m, n] = size(A);

R = max(abs(A'));    %Find the maximum element in each row
rcmax = max(R);
rcmin = min(R);

if rcmin == 0.
    fprintf('The rank of the matrix is not sufficient!\n');
    fprintf('System is not inversible!\n');
    return;
end

R = sqrt(1 ./ R);    % Invert the scale factors
eqA = A;
eqb = b;

if (rcmin/rcmax) < thresh
    fprintf('RAR scaling : yes!\n');
    eqA = spdiags(R',0,m,m)*eqA*spdiags(R',0,m,m);
    eqb = spdiags(R',0,m,m)*eqb;
else
    fprintf('RAR scaling : no\n');
    R = ones(1,m);
end

C = R;
return;

