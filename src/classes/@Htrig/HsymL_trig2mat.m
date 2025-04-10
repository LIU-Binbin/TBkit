% HSYML_TRIG2MAT Convert a symbolic HsymL_trig term to its matrix representation.
%
% SYNTAX:
%   mat = HsymL_trig2mat(H_htrig, HsymL_trig)
%
% DESCRIPTION:
%   This function converts a symbolic term (HsymL_trig) from the Htrig object into a matrix.
%   It calculates the total number of discretized slab points (NS) from the Nslab property and
%   then uses the Delta_Oper function along with coeff_extract to form a delta operator. The
%   delta operator is evaluated over all combinations of indices generated from the discretization,
%   yielding a matrix representation of the symbolic term.
%
% INPUTS:
%   H_htrig    - A Htrig object, which contains the properties Nslab and relevant methods.
%   HsymL_trig - A symbolic expression corresponding to one term in HsymL_trig.
%
% OUTPUT:
%   mat        - A double-precision matrix representation of the symbolic term.
%
% EXAMPLE:
%   % Given a Htrig object H and a symbolic term HsymL_trig_term:
%   mat = HsymL_trig2mat(H, HsymL_trig_term);
%
function mat = HsymL_trig2mat(H_htrig, HsymL_trig)
    NSLAB = (H_htrig.Nslab == 0) + H_htrig.Nslab;
    NS = prod(NSLAB);
    
    % Compute the delta operator using extracted coefficients from HsymL_trig.
    Delta_ijk = Htrig.Delta_Oper(Htrig.coeff_extract(HsymL_trig));
    
    % Generate a grid of indices over the discretized slab.
    tmpmat = meshgrid((1:NS).', (1:NS).');
    NS1 = reshape(tmpmat, NS^2, 1);
    NS2 = reshape(tmpmat.', NS^2, 1);
    
    % Convert linear indices to subscript indices.
    [i1, i2, i3] = ind2sub(NSLAB, NS1);
    [j1, j2, j3] = ind2sub(NSLAB, NS2);
    
    % Evaluate the delta operator over the index pairs and reshape the result into a matrix.
    mat = double(reshape(Delta_ijk(i1, i2, i3, j1, j2, j3), NS, NS)).';
end

