function Poly_priciplayer_mat = Poly_priciplayer_mat_gen(principle_layer)
Poly_priciplayer_mat(:,:,principle_layer) = eye(principle_layer);
for i = 1:principle_layer-1
base_mat = zeros(principle_layer);
base_mat(1:principle_layer-i,i+1:principle_layer) = eye(principle_layer-i);
Poly_priciplayer_mat(:,:,principle_layer+i) = base_mat;
Poly_priciplayer_mat(:,:,principle_layer-i) = base_mat';
end
end
