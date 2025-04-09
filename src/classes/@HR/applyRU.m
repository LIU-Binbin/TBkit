function [H_hr_R,H_hr] = applyRU(H_hr,SymOper)
%APPLYRU Apply symmetry operation to HR object (tight-binding Hamiltonian)
%   This function implements symmetry operations on tight-binding Hamiltonians
%   including handling of vector hopping terms and complex conjugation
%
%   Inputs:
%       H_hr    - HR object containing tight-binding model parameters
%       SymOper - Symmetry operation object with fields:
%                 .conjugate     : Boolean for complex conjugation
%                 .antisymmetry  : Boolean for antisymmetry operation
% 
%   Outputs:
%       H_hr_R  - Transformed HR object after symmetry operation
%       H_hr    - Original HR object (dualized version)
%
%   Throws:
%       error   - If dimension mismatch in VectorDistMat occurs
%
%   Processing flow:
%       1. Dualize input Hamiltonian using symmetry operations
%       2. Verify transformation matrix dimensions
%       3. Handle vector hopping terms with symmetry constraints
%       4. Process coefficient/numeric matrices with symmetry operations

arguments
    H_hr HR;
    SymOper     ;
end
% Dualize operator and get transformation matrix
[H_hr,VectorDistMat] = dualizeOper(H_hr,SymOper);

% Validate transformation matrix dimensions
if size(VectorDistMat,1)~=size(VectorDistMat,2)
    error('Dimension mismatch in VectorDistMat - non-square transformation matrix');
end

H_hr_R = H_hr; % Initialize output object

% Vector hopping term processing
if H_hr.vectorhopping
    % Extract vector components
    AvectorLtmp = H_hr_R.AvectorL;
    BvectorLtmp = H_hr_R.BvectorL;
    CvectorLtmp = H_hr_R.CvectorL;
    
    % Apply complex conjugation transformations
    if SymOper.conjugate
        BvectorLtmp = -(BvectorLtmp);
        CvectorLtmp(end/2+1:end,:) = -CvectorLtmp(end/2+1:end,:);
    end
    
    % Apply antisymmetry transformations
    if SymOper.antisymmetry
        AvectorLtmp = -AvectorLtmp;
        BvectorLtmp = -BvectorLtmp;
        CvectorLtmp = -CvectorLtmp;
    end
    
    % Transform components using symmetry operation
    H_hr_R.AvectorL = VectorDistMat*AvectorLtmp;
    H_hr_R.BvectorL = VectorDistMat*BvectorLtmp;
    
    % Split and transform complex components
    CL1 = VectorDistMat*CvectorLtmp(1:end/2,:);  % Real part
    CL2 = VectorDistMat*CvectorLtmp(end/2+1:end,:);  % Imaginary part
    H_hr_R.CvectorL = [real(CL1)-imag(CL2); imag(CL1)+real(CL2)];
    return;
end

% Coefficient matrix processing
if H_hr.coe
    HcoeLtmp = H_hr_R.HcoeL;
    if SymOper.conjugate
        HcoeLtmp = conj(HcoeLtmp);
    end
    if SymOper.antisymmetry
        HcoeLtmp = -HcoeLtmp;
    end
    H_hr_R.HcoeL = VectorDistMat*HcoeLtmp;
end

% Numeric matrix processing
if H_hr.num == true
    HnumLtmp = H_hr_R.HnumL;
    if SymOper.conjugate
        HnumLtmp = conj(HnumLtmp);
    end
    if SymOper.antisymmetry
        HnumLtmp = -HnumLtmp;
    end
    H_hr_R.HnumL = VectorDistMat*HnumLtmp;
end
end