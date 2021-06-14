function y = applyAI_SPD(D, N, b, e, x)

y = x;
for i = 1 : N
    y(b(i) : e(i), :) = D{i} * x(b(i) : e(i), :);
end