classdef Y_l__m < HollowKnight
properties
l   double;
m   double;
end
properties(Dependent)
expression ;
parity ;
end
properties(Dependent,Hidden)
hollow;
end
properties(Hidden)
n   = 0         ;
symbolic = true ;
explicit = false;
end
methods
 YlmObj = Y_l__m(l,m,coe,n,TH,PHI,options)
end
methods
 hollow = get.hollow(YlmObj)
 parity = get.parity(YlmObj)
 expression = get.expression(YlmObj)
end
methods
 C = plus(A,B)
 A = umius(A)
 C = minus(A,B)
 C = innertimes(A,B)
 C = eq(A,B,options)
 C = mrdivide(A,B)
 C = lt(A,B,options)
 sum()
end
methods
 B = contractrow(A,options)
end
methods
 Ak = rotateinner(A,abc,RightorLeft,immproper,conjugate,antisymmetry)
end
methods
 disp(YlmObj,options)
 expression = pretty(YlmObj,options)
 Tesseralexpansion = Tesseral(Y_lmObj,options)
 str = string(YlmObj,options)
 expression = formula(YlmObj,options)
 expression = explicitformula(YlmObj,options)
end
methods(Static)
 SymExpr = d(j,m1,m2,seed)
 SymExpr = d_mm__j(j,m1,m2,seed)
 PlmExpr = Plm(l,m,seed,options)
 NlmExpr = Nlm(l,m,seed,options)
 YlmExpr = InnerY_l__m(l,m,seed1,seed2,options)
end
methods(Static)
 DFact = DFactorial(n)
 c = binomial(n, k)
 cg = CG(j1,m1,j2,m2,j,m)
  W = w3j(j1,j2,j3,m1,m2,m3)
  W = w6j(a,b,c,d,e,f)
  W = w9j(a,b,c,d,e,f,g,h,j)
 cg = ClebschGordan(j1,j2,j,m1,m2,m)
  W = Wigner3j( j123, m123 )
  W = Wigner3j_sym(j123, m123)
end
methods(Static)
 OutExpr = nlm2atomic(l,m,n,options)
 Outstr = l__m2str(l,m)
 Outstr = lm2str(l,m)
 [dupNames, dupNdxs] = getDuplicates(aList)
 Ind = IndexOfMultiples(A)
 T = isMultiple(A)
end
methods(Static)
 tri = del(a,b,c)
 tf = triangular_cond(a,b,c)
 sym_fact = sym_fact(n)
end
end
