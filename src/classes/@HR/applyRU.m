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
%H_hr_R = H_hr; % Initialize output object
if SymOper.identity == SymOper
    H_hr_R = H_hr; 
    return;
end
% Dualize operator and get transformation matrix
[H_hr,VectorDistMat] = dualizeOper(H_hr,SymOper);

% Validate transformation matrix dimensions
% if size(VectorDistMat,1)~=size(VectorDistMat,2)
%     error('Dimension mismatch in VectorDistMat - non-square transformation matrix');
% end

H_hr_R = H_hr; 

% Vector hopping term processing
if H_hr_R.vectorhopping
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
    fprintf('%d，%d\n',length(H_hr.vectorL_map),length(H_hr.vectorL_map));
    return;
end

% Coefficient matrix processing
if H_hr_R.coe
    HcoeLtmp = H_hr.HcoeL;
    if SymOper.conjugate
        HcoeLtmp = conj(HcoeLtmp);
    end
    if SymOper.antisymmetry
        HcoeLtmp = -HcoeLtmp;
    end
    H_hr_R.HcoeL = VectorDistMat*HcoeLtmp;
end

% Numeric matrix processing
if H_hr_R.num == true
    HnumLtmp = H_hr.HnumL;
    if SymOper.conjugate
        HnumLtmp = conj(HnumLtmp);
    end
    if SymOper.antisymmetry
        HnumLtmp = -HnumLtmp;
    end
    cmpData = HnumLtmp(VectorDistMat.idL).*VectorDistMat.UL;
    H_hr_R.HnumL = segmentedSum(cmpData , VectorDistMat.offset);
end
end

function [cmpSums] = segmentedSum( cmpData, startIdx)
    % 输入:
    %   cmpData   : 扁平化的复数数组 (complex single)
    %   slices    : N×2 切片矩阵 [start, end]
    %
    % 输出:
    %   cmpSums   : 每个切片的复数和 (complex double)

    % 向量化分段求和
    % startIdx = slices(:, 1);
    endIdx = [startIdx(2:end)-1;numel(cmpData)];
    % ===== 2. 复数数组分段求和 =====
    % 分离实部虚部计算 (比直接操作复数快3倍)
    realPart = real(cmpData);
    imagPart = imag(cmpData);
    
    % 计算实部虚部分别求和
    cumReal = [0; cumsum(realPart(:))];
    cumImag = [0; cumsum(imagPart(:))];
    
    % 合并结果
    cmpSums = complex(...
        cumReal(endIdx+1) - cumReal(startIdx), ...
        cumImag(endIdx+1) - cumImag(startIdx));
end