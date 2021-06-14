function y = applyS_SPD(DNp1, RD, L, N, x)

y = DNp1 * x;
y = y - applyBtDm1B_SPD(RD, L, N, x);