function tf = triangular_cond(a,b,c)
tf = (c >= abs(a-b)) & (c <= a+b);
end
