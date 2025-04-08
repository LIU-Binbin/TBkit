function H_hk = setup_single(H_hk,symbolic_polynomial,i,j,silence)
if nargin <5
silence = false;
end
tempmat = zeros( H_hk.Basis_num);
tempmat(i,j) =1 ;
H_hk = H_hk.setup_rough(symbolic_polynomial,tempmat,silence);
end
