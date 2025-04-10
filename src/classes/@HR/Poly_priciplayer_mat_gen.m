function Poly_priciplayer_mat = Poly_priciplayer_mat_gen(principle_layer)
% POLY_PRINCIPLAYER_MAT_GEN Generate principal layer polynomial matrices
%
%   POLY_PRINCIPLAYER_MAT = POLY_PRINCIPLAYER_MAT_GEN(PRINCIPLE_LAYER)
%   generates polynomial matrices for principal layer calculations
%
%   Input:
%       principle_layer - Number of layers to generate
%   Output:
%       Poly_priciplayer_mat - Generated polynomial matrices
%
%   Notes:
%       - Creates identity matrix for principal layer
%       - Generates shifted matrices for neighboring layers
%       - Used in transport calculations
Poly_priciplayer_mat(:,:,principle_layer) = eye(principle_layer);
for i = 1:principle_layer-1
base_mat = zeros(principle_layer);
base_mat(1:principle_layer-i,i+1:principle_layer) = eye(principle_layer-i);
Poly_priciplayer_mat(:,:,principle_layer+i) = base_mat;
Poly_priciplayer_mat(:,:,principle_layer-i) = base_mat';
end
end
