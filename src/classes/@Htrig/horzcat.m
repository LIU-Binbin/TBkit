% HORZCAT Overloads the horizontal concatenation operator ( [A, B] ) for Htrig objects.
%
% SYNTAX:
%   C = horzcat(A, B)
%
% DESCRIPTION:
%   This function horizontally concatenates two Htrig objects by combining their Basis 
%   properties and forming block-diagonal matrices for the HnumL and HcoeL fields. 
%   The resultant Htrig object has a new Basis that is the vertical concatenation of the 
%   two input bases, and its Basis_num is the sum of the individual Basis_num values. 
%   For each 'kind', the corresponding HnumL and HcoeL matrices are updated using blkdiag.
%
% INPUT:
%   A, B - Htrig objects to be concatenated.
%
% OUTPUT:
%   C    - The resulting Htrig object after concatenation.
%
% EXAMPLE:
%   % Given two Htrig objects H1 and H2:
%   C = [H1, H2];
%
function C = horzcat(A, B)
    if isa(A, 'Htrig') && isa(B, 'Htrig')
        H_htrig1 = A;
        H_htrig2 = B;
        H_htrig = A;
        H_htrig.Basis = [H_htrig1.Basis; H_htrig2.Basis];
        H_htrig.Basis_num = H_htrig1.Basis_num + H_htrig2.Basis_num;
        H_htrig.HnumL = zeros(H_htrig.Basis_num, H_htrig.Basis_num, H_htrig.Kinds);
        H_htrig.HcoeL = sym(H_htrig.HnumL);
        for i = 1:H_htrig.Kinds
            H_htrig.HnumL(:, :, i) = blkdiag(H_htrig1.HnumL(:, :, i), ...
                                              H_htrig2.HnumL(:, :, i));
            H_htrig.HcoeL(:, :, i) = blkdiag(H_htrig1.HcoeL(:, :, i), ...
                                              H_htrig2.HcoeL(:, :, i));
        end
        H_htrig.Trig_to_save = sym(zeros(H_htrig.Basis_num, H_htrig.Basis_num));
    else
        % If inputs are not Htrig objects, use MATLAB's default behavior.
    end
    C = H_htrig;
end
