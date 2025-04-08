function Degree = checkDegree(Hsym,Dim)
arguments
Hsym
Dim = 3;
end
VarsSeqLcart = [sym('k_x'),sym('k_y'),sym('k_z'),sym('k_w')];
Degree_mat = polynomialDegree(Hsym,VarsSeqLcart(1:Dim));
Degree = max(Degree_mat,[],'all');
end
