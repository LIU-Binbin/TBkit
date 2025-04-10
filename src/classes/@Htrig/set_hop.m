function H_htrig = set_hop(H_htrig, SymHopping, SymVar, Indij)
%SET_HOP Add or update hopping terms in an Htrig object.
%
%   H_htrig = SET_HOP(H_htrig, SymHopping, SymVar)
%   Adds hopping matrices (SymHopping) corresponding to symbolic variables
%   (SymVar) to the Htrig object H_htrig. This form supports batch updates
%   to the Hamiltonian using either symbolic or numeric mode.
%
%   H_htrig = SET_HOP(H_htrig, SymHopping, SymVar, Indij)
%   Adds hopping terms in element (i,j) of the hopping matrix for each
%   symbolic component extracted from SymVar, scaled by the coefficients
%   in the trigonometric expression.
%
%   Inputs:
%       H_htrig     - Htrig object
%       SymHopping  - Hopping matrix/matrices (can be symbolic or numeric)
%       SymVar      - Corresponding symbolic momentum-dependent term(s)
%       Indij       - (Optional) [i, j] indices for element-wise update
%
%   Output:
%       H_htrig     - Updated Htrig object
%
%   Example:
%       H = set_hop(H, t*eye(2), [0 0 0]);           % Add onsite term
%       H = set_hop(H, t*[0 1; 1 0], [1 0 0]);       % Add hopping with phase
%       H = set_hop(H, t, sin(kx), [1 2]);           % Add trig. hopping term
%
%   Note:
%       - If Indij is provided, the function treats the term as a scalar
%         and expands trigonometric symbolic terms via split_sym_eq.
%       - The function automatically determines whether numeric or
%         symbolic mode is used based on H_htrig settings.
%
%   See also: k_symbol2Kind, add_empty_one, split_sym_eq, generalcontractrow

if isempty(H_htrig.num) || isempty(H_htrig.coe)
    [~,~,H_htrig] = H_htrig.NumOrCoe();
end

if nargin < 4
    nSymVar = size(SymVar,1);
    if size(SymVar,1) > 1 && size(SymVar,2) == 1
        SymVar = SymVar.';
    end

    if isempty(H_htrig.num) || isempty(H_htrig.coe)
        [~,~,H_htrig] = H_htrig.NumOrCoe;
    end

    if H_htrig.num
        SymHopping = double(SymHopping);
    else
        SymHopping = sym(SymHopping);
    end

    if isvector(SymHopping)
        matmode = false;
    elseif size(SymHopping,3) == 1 && ismatrix(SymHopping)
        SymHopping = repmat(SymHopping, [1 1 nSymVar]);
        matmode = true;
    else
        matmode = true;
    end

    if matmode
        for i = 1:nSymVar
            Kind = H_htrig.k_symbol2Kind(SymVar(i,:));
            if isempty(Kind)
                H_htrig = add_empty_one(H_htrig, SymVar(i));
                Kind = H_htrig.Kinds;
            else
                try
                    if Kind == 0
                        H_htrig = add_empty_one(H_htrig, SymVar(i));
                        Kind = H_htrig.Kinds;
                    end
                catch
                end
            end

            if H_htrig.coe
                H_htrig.HcoeL(:,:,Kind) = H_htrig.HcoeL(:,:,Kind) + SymHopping(:,:,i);
            else
                H_htrig.HnumL(:,:,Kind) = H_htrig.HnumL(:,:,Kind) + SymHopping(:,:,i);
            end
        end
    else
        if H_htrig.coe
            HsymL_coeLtmp = [H_htrig.HsymL_coeL; SymVar];
            Hcoetmp = [H_htrig.HcoeL; SymHopping];
            [H_htrig.HsymL_coeL, H_htrig.HcoeL] = ...
                HollowKnight.generalcontractrow(HsymL_coeLtmp, Hcoetmp);
        else
            HsymL_numLtmp = [H_htrig.HsymL_numL; SymVar];
            Hnumtmp = [H_htrig.HnumL; SymHopping];
            [H_htrig.HsymL_numL, H_htrig.HnumL] = ...
                HollowKnight.generalcontractrow(HsymL_numLtmp, Hnumtmp);
        end
    end
else
    [coeff_trig, symvar_list_trig, H_htrig] = split_sym_eq(H_htrig, simplify(SymVar));
    nc = length(coeff_trig);
    for i = 1:nc
        SymVar = symvar_list_trig(i);
        Kind = H_htrig.k_symbol2Kind(SymVar);
        if isempty(Kind)
            H_htrig = add_empty_one(H_htrig, SymVar);
            Kind = H_htrig.Kinds;
        end
        H_htrig.HcoeL(Indij(1), Indij(2), Kind) = ...
            H_htrig.HcoeL(Indij(1), Indij(2), Kind) + SymHopping * coeff_trig(i);
    end
end
end

