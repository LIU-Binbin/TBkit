function [H_hr_R, H_hr] = applyR(H_hr, R)
%APPLYR Apply spatial symmetry operation to Hamiltonian representation
%   This function implements coordinate system transformation and symmetry
%   constraint enforcement for tight-binding Hamiltonian objects
%
%   Inputs:
%       H_hr - HR object containing tight-binding parameters
%       R    - Spatial symmetry operation matrix (3x3 rotation/inversion)
%
%   Outputs:
%       H_hr_R - Transformed Hamiltonian with applied symmetry
%       H_hr    - Modified original Hamiltonian (structural changes)
%
%   Processing Flow:
%       1. Convert symmetry operation to reciprocal space representation
%       2. Prepare Hamiltonian for symmetry operation
%       3. Apply dual space transformation
%       4. Reorganize Hamiltonian parameters according to symmetry

    % Validate input types
    arguments
        H_hr HR;
        R;
    end
    
    % Convert crystal to real-space symmetry operator
    Rf = Oper.Rc2Rf(inv(R), H_hr.Rm);
    
    % Prepare Hamiltonian structure
    H_hr = H_hr.rewrite();
    
    % Apply dual space transformation and get index mapping
    [H_hr, R_vector_dist_] = dualizeR(H_hr, Rf);
    
    % Initialize output Hamiltonian
    H_hr_R = H_hr;
    
    % Reindex numerical parameters
    if H_hr.num
        H_hr_R.HnumL = H_hr_R.HnumL(R_vector_dist_);
    end
    
    % Reindex symbolic parameters
    if H_hr.coe
        H_hr_R.HcoeL = H_hr_R.HcoeL(R_vector_dist_);
    end
end