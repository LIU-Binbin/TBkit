function H_hr = set_overlap_mat(H_hr,amp,vector,mode)
V = H_hr.vectorL_overlap;
switch H_hr.Type
case 'list'
sizeamp = size(amp);
for n = 1:numel(amp)
[si,sj] = ind2sub(sizeamp,n);
if ~isequal(sym(amp(n)),sym(0))
H_hr = set_overlap_single(H_hr,amp(n),si,sj,vector,mode);
end
end
otherwise
[~,seq]=ismember((vector),V,'rows');
if seq == 0
seq = H_hr.NRPTS +1;
H_hr = H_hr.add_empty_one(vector);
end
switch mode
case 'set'
zeros_num_mat = zeros(H_hr.WAN_NUM);
if H_hr.HnumL(:,:,seq) ~= zeros_num_mat
warning('May be you should use add mode on this NRPT');
fprintf('%d %d %d %d %d \n',vector(1),vector(2),vector(3));
end
H_hr.SnumL(:,:,seq) = amp ;
case 'add'
H_hr.SnumL(:,:,seq) = H_hr.SnumL(:,:,seq) + amp  ;
case 'sym'
zeros_coe_mat = sym(zeros(H_hr.WAN_NUM));
if H_hr.ScoeL(:,:,seq) ~= zeros_coe_mat
warning('May be you should use symadd mode on this NRPT ');
fprintf('%d %d %d %d %d \n',vector(1),vector(2),vector(3));
end
H_hr.ScoeL(:,:,seq) =  amp  ;
case 'symadd'
H_hr.ScoeL(:,:,seq) = H_hr.ScoeL(:,:,seq) + amp  ;
end
end
end
