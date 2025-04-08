function C = mtimes(A,B)
if isa(A,'Htrig') && isa(B,'Htrig')
error('not be implmented');
elseif isa(A,'Htrig') && ~isa(B,'Htrig')
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
elseif ~isa(A,'Htrig') && isa(B,'Htrig')
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
