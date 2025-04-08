function J_mat = spin_matrices_from_orb(quantumL,include_0,strict)
arguments
quantumL double =[1,0,0,0.5,0.5;1,0,0,0.5,-0.5];
include_0 logical = false;
strict logical = true;
end
Norb = size(quantumL,1);
Sz =  diag(quantumL(:,5));
S_plus = zeros(Norb);
S_minus = zeros(Norb);
for i  = 1:Norb
for j = 1:Norb
J_l = quantumL(i,2)+quantumL(i,4);
jz_l = quantumL(i,5);
J_r = quantumL(j,2)+quantumL(j,4);
jz_r = quantumL(j,5);
if strict
S_plus(i,j) =  Oper.braket([J_l,jz_l ],'S+',[J_r,jz_r ]);
S_minus(i,j) =  Oper.braket([J_l,jz_l ],'S-',[J_r,jz_r ]);
else
S_plus(i,j) =  Oper.braket([J_l,jz_l ],'^S+',[J_r,jz_r ]);
S_minus(i,j) =  Oper.braket([J_l,jz_l ],'^S-',[J_r,jz_r ]);
end
end
end
Sx  = 1/2 * (S_plus+S_minus);
Sy  = 1/(2*1i) * (S_plus- S_minus);
if include_0
J_mat(:,:,1) = eye(S);
J_mat(:,:,2) = Sx;
J_mat(:,:,3) = Sy;
J_mat(:,:,4) = Sz;
else
J_mat(:,:,1) = Sx;
J_mat(:,:,2) = Sy;
J_mat(:,:,3) = Sz;
end
end
