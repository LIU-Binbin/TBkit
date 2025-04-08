function varargout = Green_prepare(H_hrz,principle_layer,fin_dir)
if nargin <3
fin_dir = 2;
end
if nargin <2
principle_layer = 3;
end
WANNUM = H_hrz.WAN_NUM;
vector_list = H_hrz.vectorL ;
label_list = vector_list(:,fin_dir);
Poly_priciplayer_mat = HR.Poly_priciplayer_mat_gen(principle_layer);
switch H_hrz.Type
case 'mat'
H00 = zeros(WANNUM*principle_layer);
for i = -(principle_layer-1):(principle_layer-1)
z = i;
label = find(label_list == z);
if label >0
H00=H00+kron(Poly_priciplayer_mat(:,:,i+principle_layer),H_hrz.HnumL(:,:,label));
end
end
H10 = zeros(WANNUM*principle_layer);
for i = 1:principle_layer
z = i;
label = find(label_list == z);
if label >0
H10=H10+kron(Poly_priciplayer_mat(:,:,i),H_hrz.HnumL(:,:,label));
end
end
case 'sparse'
H00 = sparse(WANNUM*principle_layer);
for i = -(principle_layer-1):(principle_layer-1)
z = i;
label = find(label_list == z);
if label >0
H00=H00+kron(Poly_priciplayer_mat(:,:,i+principle_layer),H_hrz.HnumL{label});
end
end
H10 = sparse(WANNUM*principle_layer);
for i = 1:principle_layer
z = i;
label = find(label_list == z);
if label >0
H10=H10+kron(Poly_priciplayer_mat(:,:,i),H_hrz.HnumL{label});
end
end
end
H01 = H10' ;
if nargout == 2
varargout{1} = H00;
varargout{2} = H01;
return;
else
varargout{1} = H00;
varargout{2} = H01;
end
switch H_hrz.Type
case 'mat'
H_hrz_green.HnumL(:,:,1) = H10;H_hrz_green.vectorL(1,:) = ([0 ,0, -1]);
H_hrz_green.HnumL(:,:,2) = H00;H_hrz_green.vectorL(2,:) = ([0 ,0, 0]);
H_hrz_green.HnumL(:,:,3) = H01;H_hrz_green.vectorL(3,:) = ([0 ,0, 1]);
case 'sparse'
H_hrz_green.HnumL{1} = H10;H_hrz_green.vectorL(1,:) = ([0 ,0, -1]);
H_hrz_green.HnumL{2}= H00;H_hrz_green.vectorL(2,:) = ([0 ,0, 0]);
H_hrz_green.HnumL{3} = H01;H_hrz_green.vectorL(3,:) = ([0 ,0, 1]);
end
varargout{3} = H_hrz_green;
end
