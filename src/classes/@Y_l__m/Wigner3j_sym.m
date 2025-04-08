function  W = Wigner3j_sym(j123, m123)
a = j123(1); b = j123(2); c = j123(3);
alf = m123(1); bet = m123(2); gam = m123(3);
tj_iszero = false;
if abs(alf) > abs(a) || abs(bet) > abs(b) || abs(gam) > abs(c)
tj_iszero = true;
end
if alf + bet + gam ~= 0
tj_iszero = true;
end
if abs(c) > abs(a + b) || abs(c) < abs(a - b)
tj_iszero = true;
end
if ~tj_iszero
tri_coef = Y_l__m.sym_fact(a + b - c)*Y_l__m.sym_fact(a - b + c)*Y_l__m.sym_fact(-a + b + c)/Y_l__m.sym_fact(a + b + c + 1);
pre_coef = Y_l__m.sym_fact(a + alf)*Y_l__m.sym_fact(a - alf)...
*Y_l__m.sym_fact(b + bet)*Y_l__m.sym_fact(b - bet)...
*Y_l__m.sym_fact(c + gam)*Y_l__m.sym_fact(c - gam);
tmax = max([a+alf,a-alf, b+bet,b-bet, c+gam,c-gam,...
a+b-c,b+c-a,c+a-b]);
t_sum = sym(0);
for t = 0:tmax
if c-b+t+alf >= 0 && c-a+t-bet >=0 && a+b-c-t >= 0 && a-t-alf >=0 && b-t+bet >= 0
t_denom = Y_l__m.sym_fact(t)*Y_l__m.sym_fact(a+b-c-t)*...
Y_l__m.sym_fact(c-b+t+alf)*Y_l__m.sym_fact(c-a+t-bet)*...
Y_l__m.sym_fact(a-t-alf)*Y_l__m.sym_fact(b-t+bet);
t_sum = t_sum + sym((-1)^t)/t_denom;
end
end
W = sym((-1)^(a-b-gam))*sqrt(tri_coef*pre_coef)*t_sum;
end
if tj_iszero
W = sym(0);
end
W = simplify(W);
end
