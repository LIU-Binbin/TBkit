function H_hk = applyU(H_hk, U, conjugate, antisymmetry)
%APPLYU Apply unitary transformation to HK object
%   H_hk = APPLYU(H_hk, U, conjugate, antisymmetry) applies a unitary
%   transformation matrix U to the HK object's coefficients, with optional
%   complex conjugation and antisymmetry operations.
%
% Input Arguments:
%   H_hk        : HK object
%       Input Hamiltonian object
%   U           : double/sym matrix
%       Unitary transformation matrix (default: nan)
%   conjugate   : logical
%       Flag for complex conjugation (default: false)
%   antisymmetry: logical
%       Flag for sign inversion (default: false)
%
% Output Arguments:
%   H_hk : HK object
%       Transformed HK object
%
% Methods Called:
%   HK.matrixtimespage    - Page-wise matrix multiplication
%   HK.page_mtimes_matrix - Matrix-page multiplication
%
% Example:
%   U = rand(3); 
%   Hk_transformed = applyU(Hk, U, true, false);

arguments
    H_hk HK
    U = nan
    conjugate logical = false
    antisymmetry logical = false
end

% Handle Oper object input
if isa(U, 'Oper')
    conjugate = U.conjugate;
    antisymmetry = U.antisymmetry;
    U = U.U;
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

% Calculate inverse matrix
U_inv = inv(U);

% Process symbolic coefficients
if coe_label
    HcoeLtmp = H_hk.HcoeL;
    if conjugate
        HcoeLtmp = conj(HcoeLtmp);
        HcoeLtmp = HK.matrixtimespage(HK.factorlist_parity(H_hk.Degree), HcoeLtmp);
    end
    if antisymmetry
        HcoeLtmp = -HcoeLtmp;
    end
    HcoeLtmp = HK.page_mtimes_matrix(HK.matrix_mtimes_page(U, HcoeLtmp), U_inv);
    H_hk.HcoeL = HcoeLtmp;
end

% Process numerical coefficients
if num_label
    U_page = repmat(U, [1 1 size(H_hk.HnumL, 3)]);
    U_inv_page = repmat(U_inv, [1 1 size(H_hk.HnumL, 3)]);
    HnumLtmp = H_hk.HnumL;
    if conjugate
        HnumLtmp = conj(HnumLtmp);
        HnumLtmp = HK.matrixtimepage(HK.factorlist_parity(H_hk.Degree), HnumLtmp);
    end
    if antisymmetry
        HnumLtmp = -HnumLtmp;
    end
    HnumLtmp = pagemtimes(pagemtimes(U_page, HnumLtmp), U_inv_page);
    H_hk.HnumL = HnumLtmp;
end
end
