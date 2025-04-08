function H_htrig = translate(H_htrig,U)
U_inv = inv(U);
if isa(U,'double')
for i = 1:C.Kinds
H_htrig.HcoeL(:,:,i) = U_inv*H_htrig.HcoeL(:,:,i)*U;
H_htrig.HnumL(:,:,i) = U_inv*H_htrig.HnumL(:,:,i)*U;
end
elseif isa(U,'sym')
for i = 1:H_htrig.Kinds
H_htrig.HcoeL(:,:,i) = U_inv*H_htrig.HcoeL(:,:,i)*U;
end
else
end
end
