function C = mtimes_inv(B,A)
if ~isa(A,'HK') && isa(B,'HK')
C = B;
if isa(A,'double')
for i = 1:C.Kinds
C.HcoeL(:,:,i) = A*C.HcoeL(:,:,i);
C.HnumL(:,:,i) = A*C.HnumL(:,:,i);
end
elseif isa(A,'sym')
for i = 1:C.Kinds
C.HcoeL(:,:,i) = A*C.HcoeL(:,:,i);
end
else
end
end
end
