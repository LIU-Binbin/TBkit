function C = mrdivide(A,B)
if isa(A,'HK') && isa(B,'HK')
elseif isa(A,'HK') && ~isa(B,'HK')
C = A;
C.HcoeL =  C.HcoeL/B;
else
end
end
