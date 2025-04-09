function H_hk = setup(H_hk, Var_cell, k_cell, mat_cell, silence)
%SETUP Configure Hamiltonian terms from cell arrays
%
% Syntax:
%   Hk = setup(Hk, Var_cell, k_cell, mat_cell)
%   Hk = setup(..., silence)
%
% Inputs:
%   Hk - HK Hamiltonian to configure
%   Var_cell - Cell array of coefficient variables
%   k_cell - Cell array of k-space symbols
%   mat_cell - Cell array of matrix elements
%   silence - Suppress progress if true (default=false)
%
% Output:
%   Hk - Configured Hamiltonian
%
% Description:
%   Core method for building Hamiltonian by:
%   1. Mapping k-space terms to matrix elements
%   2. Handling both symbolic and numeric coefficients
%   3. Supporting batch processing of multiple terms
%
% Note:
%   Uses k_symbol2Kind for term identification
%   Progress bar shown unless silenced
%
% Example:
%   Hk = Hk.setup({sym('a')}, {'k_x'}, {sigma_x});
if nargin < 5
    silence = false;
end
if length(Var_cell)~=length(k_cell) && length(k_cell)~=length(mat_cell)
    error('error!');
end
if ~silence
    pb = TBkit_tool_outer.CmdLineProgressBar('Setting ');
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
