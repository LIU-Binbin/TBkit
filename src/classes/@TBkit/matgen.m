function matcell = matgen(H_hk,R,Accuracy)
%MATGEN Generate matrix representation of Hamiltonian for given R
%   MATCELL = MATGEN(H_hk, R, Accuracy) generates matrix representation of
%   Hamiltonian terms for given lattice vector R.
%
%   Inputs:
%       H_hk    - HK object representing the Hamiltonian
%       R       - Lattice vector (can be symbolic)
%       Accuracy - Numerical rounding accuracy (default: 6)
%
%   Output:
%       matcell - Cell array containing matrix representations for each order
%
%   See also HK
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
