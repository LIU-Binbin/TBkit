function H_hk = setup(H_hk,Var_cell,k_cell,mat_cell,silence)
%SETUP Configure HK object with Hamiltonian terms
%   H_HK = SETUP(H_HK, VAR_CELL, K_CELL, MAT_CELL, SILENCE) sets up the
%   Hamiltonian terms in an HK object.
%
%   Inputs:
%       H_hk     - Input HK object
%       Var_cell - Cell array of variables/coefficients
%       k_cell   - Cell array of k-point symbols
%       mat_cell - Cell array of matrix elements
%       silence  - Suppress progress display if true (default: false)
%
%   Output:
%       H_hk     - Configured HK object
if nargin < 5
    silence = false;
end
if length(Var_cell)~=length(k_cell) && length(k_cell)~=length(mat_cell)
    error('error!');
end
if ~silence
    pb = CmdLineProgressBar('Setting ');
end
nVar_cell = length(Var_cell);
for i =1:nVar_cell
    Var = Var_cell{i};
    k_symbol = k_cell{i};
    if isa(mat_cell{i},'sym')
        matcell = sum(sym(mat_cell{i}),3);
    else
        matcell = sum(double(mat_cell{i}),3);
    end
    Kind = H_hk.k_symbol2Kind(k_symbol);
    if ~silence
        pb.print(i,nVar_cell,' term into HK ...');
        %                      fprintf('setting %dth (%s) Hk\n',Kind,string(H_hk.HsymL(Kind)));
    end
    switch class(Var)
        case 'sym'
            H_hk.HcoeL(:,:,Kind) = H_hk.HcoeL(:,:,Kind)+ matcell*Var;
        case 'string'
        case 'double'
            H_hk.HnumL(:,:,Kind) = H_hk.HcoeL(:,:,Kind)+ matcell*Var;
    end
end
if ~silence
    pb.delete();
end
end
