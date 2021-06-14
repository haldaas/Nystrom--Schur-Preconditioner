function y = applyU_SPD(L, b, e, N, x)

% determin the length of x
s = 0;
for i = 1 : N
    s = s + size(L{i}, 2);
end

y = zeros(s, size(x, 2));
for i = 1 : N
    y(b(i) : e(i), :) = L{i}' * x;
end