function [H_htrig,H_htrig2] = applyR(H_htrig,R)
% APPLYR Applies transformation matrix R to Htrig object components.
%   [H_htrig, H_htrig2] = APPLYR(H_htrig, R) transforms the numerical/coefficient
%   matrices of the Htrig object using transformation matrix R, preserving the
%   original object in H_htrig2. Operates on both numerical and symbolic coefficient
%   components based on their availability.
%
%   Inputs:
%       H_htrig - Target Htrig object to transform
%       R       - Transformation matrix (typically 3x3 spatial transformation)
%
%   Outputs:
%       H_htrig  - Transformed Htrig object with updated HnumL/HcoeL matrices
%       H_htrig2 - Original Htrig object copy before transformation
%
%   Key Operations:
%       1. Generates transformation matrix Smat via Smatgen
%       2. Applies Smat to numerical terms (HnumL) if available
%       3. Applies Smat to symbolic coefficients (HcoeL) if available
%       4. Preserves original object state in second output
%
%   Example:
%       % Apply 90-degree rotation to tight-binding model
%       R = [0 -1 0; 1 0 0; 0 0 1]; % z-axis rotation
%       [H_transformed, H_original] = applyR(H_htrig, R);
%
%   See also: Smatgen, matrixtimespage, applyOper

% Type validation
arguments
    H_htrig Htrig;
    R ;
end

% Determine numerical/symbolic component status
[num_label,coe_label] = H_htrig.NumOrCoe();

% Generate transformation matrix and preserve original state
[H_htrig,Smat] = H_htrig.Smatgen(R);
H_htrig2 = H_htrig; % Preserve pre-transformation state

% Apply transformation to numerical components
if num_label
    H_htrig.HnumL = matrixtimespage(Smat,H_htrig.HnumL);
end

% Apply transformation to symbolic coefficients
if coe_label
    H_htrig.HcoeL = matrixtimespage(Smat,H_htrig.HcoeL);
end
end
