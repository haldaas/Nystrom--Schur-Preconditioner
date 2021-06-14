function y = applyL_SPD(L, b, e, N, x)

y = zeros(size(L{1}, 1), size(x, 2));
for i = 1 : N
    y = y + L{i} * x(b(i) : e(i), :);
end