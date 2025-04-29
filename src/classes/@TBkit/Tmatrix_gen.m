function Tmatrix = Tmatrix_gen(H00,H01,w,eta)
%TMATRIX_GEN Compute transfer matrix for Green's function calculation
%
%   Syntax:
%       Tmatrix = Tmatrix_gen(H00,H01,w,eta)
%
%   Description:
%       Constructs transfer matrix for recursive Green's function method.
%       Handles ill-conditioned cases with numerical stabilization.
%
%   Inputs:
%       H00 - Intra-layer Hamiltonian
%       H01 - Inter-layer coupling
%       w   - Energy/frequency point
%       eta - Broadening parameter
%
%   Output:
%       Tmatrix - Computed transfer matrix
%
%   See also: GW_iter, GREENCAR_gen
Dimi = length(H00);
wc = w + 1i*eta;
W = (wc*eye(Dimi)-H00);

if rcond(H01) < 1e-10

    %temp_diag = diag(H01)==0;
    %disp('work');

    H01 = H01 +  eta*eye(Dimi);
    %isnan(H01)

    %H01 = (10^-nfin)*eye(Dimi);
    %disp('nor work');
    %Tmatrix = [pinv(H01)*W,-pinv(H01)*(H01');eye(Dimi) zeros(Dimi)];
    Tmatrix = [H01\W,-H01\(H01');eye(Dimi) zeros(Dimi)];
    if rcond(Tmatrix) < 1e-10
        Tmatrix = [lsqminnorm(H01,W),-lsqminnorm(H01,(H01'));eye(Dimi) zeros(Dimi)];
        return;
    end
else
    Tmatrix = [H01\W,-H01\(H01');eye(Dimi) zeros(Dimi)];
end


end