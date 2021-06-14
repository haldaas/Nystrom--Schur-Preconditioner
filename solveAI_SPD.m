function y = solveAI_SPD(RD, N, b, e, x)

for i = 1 : N
    y(b(i):e(i), :) = RD{i}\(RD{i}'\x(b(i) : e(i), :));
end