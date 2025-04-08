function H_hr = add_empty_one(H_hr,vector)
if H_hr.vectorhopping
H_hr.vectorL = [H_hr.vectorL;vector];
nvector  = size(vector,1);
H_hr.AvectorL = blkdiag(H_hr.AvectorL ,eye(nvector));
H_hr.BvectorL = blkdiag(H_hr.BvectorL ,eye(nvector));
H_hr.CvectorL = [blkdiag(H_hr.CvectorL(1:end/2,:),1*eye(nvector));...
blkdiag(H_hr.CvectorL(end/2+1:end,:),1*eye(nvector))];
return;
end
for i = 1:size(vector,1)
vector_single = (vector(i,:));
try
if (ismember(vector_single,H_hr.vectorL,'rows') && ~H_hr.overlap) || ...
(ismember(vector_single,H_hr.vectorL,'rows') ...
&& ismember(vector_single,H_hr.vectorL_overlap,'rows') && H_hr.overlap)
continue;
end
catch
end
NRPTS_new = H_hr.NRPTS +1;
H_hr.vectorL(NRPTS_new,:) = (vector_single);
if strcmp(H_hr.Type ,'mat')
if H_hr.coe
H_hr.HcoeL(:,:,NRPTS_new) = sym(zeros(H_hr.WAN_NUM));
end
if H_hr.num
H_hr.HnumL(:,:,NRPTS_new) = zeros(H_hr.WAN_NUM);
end
if H_hr.overlap
H_hr.vectorL_overlap(NRPTS_new,:) = vector_single;
if H_hr.coe
H_hr.ScoeL(:,:,NRPTS_new) = sym(zeros(H_hr.WAN_NUM));
end
if H_hr.num
H_hr.SnumL(:,:,NRPTS_new) = zeros(H_hr.WAN_NUM);
end
end
elseif strcmp(H_hr.Type ,'sparse')
if H_hr.coe
H_hr.HcoeL{NRPTS_new} = sym(sparse(H_hr.WAN_NUM));
end
if H_hr.num
H_hr.HnumL{NRPTS_new} = sparse(H_hr.WAN_NUM);
end
if H_hr.overlap
H_hr.vectorL_overlap(NRPTS_new,:) = vector_single;
if H_hr.coe
H_hr.ScoeL{NRPTS_new} = sym(sparse(H_hr.WAN_NUM));
end
if H_hr.num
H_hr.SnumL{NRPTS_new} = sparse(H_hr.WAN_NUM);
end
end
elseif strcmp(H_hr.Type ,'list')
if  H_hr.coe
H_hr.HcoeL(NRPTS_new,1) = sym(0);
end
if H_hr.num
H_hr.HnumL(NRPTS_new,1) = 0;
end
if H_hr.overlap
H_hr.ScoeL(NRPTS_new,1) = sym(0);
H_hr.SnumL(NRPTS_new,1) = 0;
end
end
end
end
