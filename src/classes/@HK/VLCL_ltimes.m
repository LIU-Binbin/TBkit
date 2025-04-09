function [vector_list, Coeffs_list] = VLCL_ltimes(VL1, CL1, VL2, CL2)
%VLCL_LTIMES Combine vector and coefficient lists
%
% Syntax:
%   [v_out, c_out] = VLCL_ltimes(VL1, CL1, VL2, CL2)
%
% Inputs:
%   VL1 - First vector list (N×3)
%   CL1 - First coefficient list (N×1)
%   VL2 - Second vector list (M×3)
%   CL2 - Second coefficient list (M×1)
%
% Outputs:
%   vector_list - Combined vectors (N*M×3)
%   Coeffs_list - Combined coefficients (N*M×1)
%
% Description:
%   Computes all combinations of input vectors/coefficients:
%   - Output vectors are element-wise sums
%   - Output coefficients are products
%
% Note:
%   Core utility for tight-binding mapping
%   Preserves relative ordering of combinations
%
% Example:
%   [v,c] = VLCL_ltimes([1 0 0],2,[0 1 0],3) % Returns [1 1 0],6
vector_list = zeros(size(CL1,1)*size(CL2,1),3);
Coeffs_list = zeros(size(CL1,1)*size(CL2,1),1);
count = 0;
for i = 1:size(VL1,1)
    for j = 1:size(VL2,1)
        count = count +1;
        vector_list(count,:) = VL1(i,:)+VL2(j,:);
        Coeffs_list(count) = CL1(i)*CL2(j);
    end
end
end
