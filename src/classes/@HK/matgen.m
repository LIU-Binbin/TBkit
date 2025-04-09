function matcell = matgen(H_hk,R,Accuracy)
%MATGEN Generate transformation matrices for position-space operations
%
% Syntax:
%   matcell = matgen(H_hk,R)
%   matcell = matgen(H_hk,R,Accuracy)
%
% Inputs:
%   H_hk - HK Hamiltonian object
%   R - Position vector or transformation matrix
%   Accuracy - Rounding precision for numeric output (default=6)
%
% Output:
%   matcell - Cell array of transformation matrices for each order
%
% Description:
%   Generates matrices that transform k-space terms under position-space
%   operation R by:
%   1. Substituting k → R·k in each term
%   2. Re-expanding in original basis
%   3. Extracting transformation coefficients
%
% Note:
%   Automatically handles both symbolic and numeric R inputs
%   Uses HK.orderlist to maintain term ordering
%
% Example:
%   mats = Hk.matgen([1 0 0; 0 -1 0; 0 0 -1]); % Mirror transformation
arguments
    H_hk HK;
    R  ;
    Accuracy double = 6;
end
if isa(R,'sym')
    coe_label = true;
else
    coe_label = false;
end
VarUsing = H_hk.VarsSeqLcart(1:H_hk.Dim);
for i = 1:H_hk.Degree+1
    Orderlist = HK.orderlist(i-1);
    nOrderlist = length(Orderlist);
    if coe_label
        matcell{i} = sym(zeros(nOrderlist));
    else
        matcell{i} = zeros(nOrderlist);
    end
    HsymC{i}(1:nOrderlist) = H_hk.HsymL_k(Orderlist);
end
HsymC_bk = HsymC;
k_orgin = VarUsing.';
k_R = R*k_orgin;
for i = 1:H_hk.Degree+1
    HsymC{i} = subs(HsymC{i},k_orgin,k_R);
end
for i = 1:H_hk.Degree+1
    H_symL_i = simplify(HsymC{i}) ;
    nOrderlist = length(H_symL_i);
    for j = 1:nOrderlist
        [A,B] = coeffs(H_symL_i(j),VarUsing);
        for k = 1:length(B)
            tempSym = B(k);
            for l  = 1:nOrderlist
                if isequal(tempSym,HsymC_bk{i}(l))
                    matcell{i}(j,l)=A(k);
                    break;
                end
            end
        end
    end
    if ~coe_label
        matcell{i} = roundn(matcell{i},-Accuracy);
    end
end
end
