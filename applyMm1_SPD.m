function y = applyMm1_SPD(L, Dm1, tildeSm1, pb, pe, N, y1, y2)
t1 = y1;
t2 = y2 - applyL_SPD(L, pb, pe, N, Dm1(y1));

t1 = Dm1(t1);
t2 = tildeSm1(t2);

y = [t1 - Dm1(applyU_SPD(L, pb, pe, N, t2)); t2];
