function S_mat = spin_matrices(s, include_0)
arguments
s double {Oper.mustBeHalfInteger(s)} =1/2;
include_0 logical = false;
end
S = 2*s+1;
Sz= 1/2 * diag(S-1:-2:-S);
S_plus = zeros(S-1,1);
for i  = 1:S-1
S_plus(i) =  sqrt((s+1)*2*i-i*(i+1))/2 ;
end
Sx  = diag(S_plus,1) +diag(S_plus,-1);
Sy  = -1i*diag(S_plus,1) +1i*diag(S_plus,-1);
if include_0
S_mat(:,:,1) = eye(S);
S_mat(:,:,2) = Sx;
S_mat(:,:,3) = Sy;
S_mat(:,:,4) = Sz;
else
S_mat(:,:,1) = Sx;
S_mat(:,:,2) = Sy;
S_mat(:,:,3) = Sz;
end
end
