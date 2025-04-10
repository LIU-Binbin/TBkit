% DUALIZE Compute the dual representation for a Htrig object (for 'exp' type).
%
%   H_htrig = dualize(H_htrig) examines the conjugate of each element in 
%   HsymL_trig for Htrig objects with Type 'exp'. For each conjugated symbol,
%   it checks its existence using k_symbol2Kind. If the symbol is not already 
%   present, it is added as a new "kind" and the corresponding coefficient (HcoeL)
%   and numerical (HnumL) matrices are initialized to zero. Finally, the 
%   Duality_vector_dist property is updated to indicate the positions of the dual symbols.
%
%   Input:
%       H_htrig - A Htrig object.
%
%   Output:
%       H_htrig - The updated Htrig object with a dual representation.
%
%   Example:
%       H_htrig = dualize(H_htrig);
%

function H_htrig = dualize(H_htrig)
    if strcmp(H_htrig.Type, 'exp')
        BASIS_NUM = H_htrig.Basis_num;
        for i = 1:length(H_htrig.HsymL_trig)
            k_symbol = conj(H_htrig.HsymL_trig(i));
            Kind = H_htrig.k_symbol2Kind(k_symbol);
            if isempty(Kind)
                Kind = H_htrig.Kinds + 1;
                H_htrig.HsymL_trig(Kind) = k_symbol;
                H_htrig.HcoeL(:, :, Kind) = sym(zeros(BASIS_NUM, BASIS_NUM, 1));
                H_htrig.HnumL(:, :, Kind) = zeros(BASIS_NUM, BASIS_NUM, 1);
            end
        end
        [~, H_htrig.Duality_vector_dist] = ismember(conj(H_htrig.HsymL_trig), H_htrig.HsymL_trig);
    end
end
