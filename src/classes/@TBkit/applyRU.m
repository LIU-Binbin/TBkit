% applyRU.m
function [H_hk_R, H_hk] = applyRU(H_hk, SymOper)
%APPLYRU Apply symmetry operation to HK object with R/U transformations
%   [H_hk_R, H_hk] = APPLYRU(H_hk, SymOper) applies a combined symmetry operation 
%   defined by SymOper structure to the HK object. The operation includes:
%   - R matrix transformation (direct/inverse based on conjugate flag)
%   - U matrix transformation (if provided)
%   - Conjugation and antisymmetry operations
%
% Input Arguments:
%   H_hk    : HK object
%       Input Hamiltonian object with coefficient matrices
%   SymOper : struct
%       Symmetry operation structure containing fields:
%       - R: 3x3 transformation matrix (double/sym)
%       - conjugate: logical flag for complex conjugation
%       - antisymmetry: logical flag for sign inversion
%       - U: Unitary matrix for basis transformation (optional)
%
% Output Arguments:
%   H_hk_R : HK object
%       Transformed HK object after R/U operations
%   H_hk   : HK object
%       Original input HK object (unmodified)
%
% Methods Called:
%   HK.matgen          - Generate matrix cell array
%   HK.orderlist       - Get polynomial order indices
%   HK.matrixtimespage - Page-wise matrix multiplication
%   HK.page_mtimes_matrix - Matrix-page multiplication
%
% Example:
%   SymOp = struct('R', eye(3), 'conjugate', false, 'antisymmetry', false);
%   [Hk_transformed, Hk_orig] = applyRU(Hk, SymOp);

% Early return for identity transformation
if ~SymOper.conjugate && ~SymOper.antisymmetry && isequal(SymOper.R, eye(3))
    H_hk_R = H_hk;
    return;
end

% Check numerical coefficients processing flag
if isequal(zeros(size(H_hk.HnumL)), H_hk.HnumL)
    num_label = false;
else
    num_label = true;
end

% Check symbolic coefficients processing flag
if isequal(sym(zeros(size(H_hk.HcoeL))), H_hk.HcoeL)
    coe_label = false;
else
    coe_label = true;
end

% Calculate inverse R matrix with sign based on conjugate flag
if SymOper.conjugate
    R_inv = -inv(SymOper.R);
else
    R_inv = inv(SymOper.R);
end

% Initialize output object and generate transformation matrices
H_hk_R = H_hk;
matcell = H_hk.matgen(R_inv);

% Process numerical coefficients with R transformation
if num_label
    for i = 0:H_hk.Degree
        Orderlist = HK.orderlist(i);
        H_hk_R.HnumL(:,:,Orderlist) = HK.matrixtimespage(...
            matcell{i+1}.', H_hk_R.HnumL(:,:,Orderlist));
    end
end

% Process symbolic coefficients with R transformation
if coe_label
    for i = 0:H_hk.Degree
        Orderlist = HK.orderlist(i);
        H_hk_R.HcoeL(:,:,Orderlist) = HK.matrixtimespage(...
            matcell{i+1}.', H_hk_R.HcoeL(:,:,Orderlist));
    end
end

% Process U transformation if provided
if ~isnan(SymOper.U)
    U = SymOper.U;
    U_inv = inv(U);
    
    % Apply U transformation to symbolic coefficients
    if coe_label
        HcoeLtmp = H_hk_R.HcoeL;
        if SymOper.conjugate
            HcoeLtmp = conj(HcoeLtmp);
        end
        if SymOper.antisymmetry
            HcoeLtmp = -HcoeLtmp;
        end
        HcoeLtmp = HK.page_mtimes_matrix(HK.matrix_mtimes_page(U, HcoeLtmp), U_inv);
        H_hk_R.HcoeL = HcoeLtmp;
    end
    
    % Apply U transformation to numerical coefficients
    if num_label
        U_page = repmat(U, [1 1 size(H_hk_R.HnumL, 3)]);
        U_inv_page = repmat(U_inv, [1 1 size(H_hk_R.HnumL, 3)]);
        HnumLtmp = H_hk_R.HnumL;
        if SymOper.conjugate
            HnumLtmp = conj(HnumLtmp);
            HnumLtmp = HK.matrixtimepage(HK.factorlist_parity(H_hk_R.Degree), HnumLtmp);
        end
        if SymOper.antisymmetry
            HnumLtmp = -HnumLtmp;
        end
        HnumLtmp = pagemtimes(pagemtimes(U_page, HnumLtmp), U_inv_page);
        H_hk_R.HnumL = HnumLtmp;
    end
end
end
