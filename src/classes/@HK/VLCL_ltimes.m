function [vector_list,Coeffs_list] = VLCL_ltimes(VL1,CL1,VL2,CL2)
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
