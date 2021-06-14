function y = applyBtDm1B_SPD(RD, L, N, x)

y = zeros(size(L{1}, 1), size(x, 2));
for i = 1 : N
    t = L{i}' * x;
    y = y + L{i} * (RD{i}\(RD{i}'\t));
end