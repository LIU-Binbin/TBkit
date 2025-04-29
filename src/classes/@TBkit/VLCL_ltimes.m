function [vector_list,Coeffs_list] = VLCL_ltimes(VL1,CL1,VL2,CL2)
%VLCL_LTIMES Combine vector lists and coefficient lists
%
%   Syntax:
%       [vector_list,Coeffs_list] = VLCL_ltimes(VL1,CL1,VL2,CL2)
%
%   Description:
%       Computes all pairwise combinations of vectors and coefficients
%       from two input lists.
%
%   Inputs:
%       VL1 - First vector list (N1 x 3)
%       CL1 - First coefficient list (N1 x 1)
%       VL2 - Second vector list (N2 x 3)
%       CL2 - Second coefficient list (N2 x 1)
%
%   Outputs:
%       vector_list - Combined vector list (N1*N2 x 3)
%       Coeffs_list - Combined coefficient list (N1*N2 x 1)
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
