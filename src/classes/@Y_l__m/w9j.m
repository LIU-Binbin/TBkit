function  W = w9j(a,b,c,d,e,f,g,h,j)
assert(isscalar(a) && isscalar(b) && isscalar(c) && ...
isscalar(d) && isscalar(e) && isscalar(f) && ...
isscalar(g) && isscalar(h) && isscalar(j),'All inputs must be scalars.')
args = [a b c d e f g h j];
if any(args < 0) || any(mod(2*args,1) ~= 0)
error('All arguments to Wigner 9-j symbol must be integer or half-integer non-negative numbers.')
elseif ~triangular_cond(a,b,c) || ~triangular_cond(d,e,f) || ...
~triangular_cond(g,h,j) || ~triangular_cond(a,d,g) || ...
~triangular_cond(b,e,h) || ~triangular_cond(c,f,j)
W = 0;
return
end
xlims(1) = a+j;
xlims(2) = b+f;
xlims(3) = d+h;
xlims(4) = a-j;
xlims(5) = b-f;
xlims(6) = d-h;
W = 0;
for x = max(xlims(4:6)):min(xlims(1:3))
W = W + (-1)^(2*x)*(2*x+1)*...
w6j(a,b,c,f,j,x)*w6j(d,e,f,b,x,h)*w6j(g,h,j,x,a,d);
end
end
