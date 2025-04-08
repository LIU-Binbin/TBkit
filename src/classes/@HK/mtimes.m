function C = mtimes(A,B)
if isa(A,'HK') && isa(B,'HK')
C = A;
elseif isa(A,'HK') && ~isa(B,'HK')
C =A;
if isa(B,'double')
for i = 1:C.Kinds
C.HcoeL(:,:,i) = C.HcoeL(:,:,i)*B;
C.HnumL(:,:,i) = C.HnumL(:,:,i)*B;
end
elseif isa(B,'sym')
for i = 1:C.Kinds
C.HcoeL(:,:,i) = C.HcoeL(:,:,i)*B;
end
else
end
elseif ~isa(A,'HK') && isa(B,'HK')
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
else
C = 0;
end
end
