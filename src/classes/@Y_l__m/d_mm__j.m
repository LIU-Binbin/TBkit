function SymExpr = d_mm__j(j,m1,m2,seed)
arguments
    j double{mustBePositive};
    m1 double
    m2 double
    seed sym = sym('theta');
end
if m1 <0 &&  m2 <0
    SymExpr = Y_l__m.d_mm__j(j,-m2,-m1,seed);
    return;
end
theta = seed;
k = min([j+m2,j-m2,j+m1,j-m1]);
switch k
    case j+m2
        a = m1 - m2;lambda = m1-m2;
    case j-m2
        a = m2 - m1;lambda = 0;
    case j+m1
        a = m2 - m1;lambda = 0;
    case j-m1
        a = m1 - m2;lambda = m1-m2;
end
b = 2*j-2*k-a;
SymExpr = (-1)^lambda*...
    nchoosek(2*j-k,k+a)^0.5*nchoosek(k+b,b)^(-0.5)*...
    sin(theta/2)^a* cos(theta/2)^b*...
    jacobiP(k,a,b,cos(theta));
end
