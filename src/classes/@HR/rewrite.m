function H_hr = rewrite(H_hr,options)
arguments
H_hr HR;
options.rewind = false;
options.Accuracy = 1e-6;
options.type = '';
end
WANNUM = H_hr.WAN_NUM;
if ~strcmp(H_hr.Type ,'list')
if  H_hr.num
if isvector(H_hr.HnumL)
warning('May not need to rewrite. Do nothing');
return;
end
NRPTS_ = numel(H_hr.HnumL);
sizeHcoeL = size(H_hr.HnumL);
HnumLtmp = reshape(H_hr.HnumL,[NRPTS_,1]);
[iL,jL,nL]= ind2sub(sizeHcoeL,1:NRPTS_);
vectorList = [H_hr.vectorL(nL,:),iL.',jL.'];
H_hr.HnumL = HnumLtmp;
H_hr.vectorL = vectorList;
if H_hr.overlap
NRPTS_S = numel(H_hr.SnumL);
sizeScoeL = size(H_hr.SnumL);
SnumLtmp = reshape(H_hr.SnumL,[NRPTS_S,1]);
[iL,jL,nL]= ind2sub(sizeScoeL,1:NRPTS_S);
vectorList_overlap = [H_hr.vectorL_overlap(nL,:),iL.',jL.'];
H_hr.SnumL = SnumLtmp;
H_hr.vectorL_overlap = vectorList_overlap;
end
elseif H_hr.coe && ~H_hr.num
if isvector(H_hr.HcoeL)
warning('May not need to rewrite. Do nothing');
return;
end
NRPTS_ = numel(H_hr.HcoeL);
sizeHcoeL = size(H_hr.HcoeL);
HcoeLtmp = reshape(H_hr.HcoeL,[NRPTS_,1]);
[iL,jL,nL]= ind2sub(sizeHcoeL,1:NRPTS_);
vectorList = [H_hr.vectorL(nL,:),iL.',jL.'];
H_hr.HcoeL = HcoeLtmp;
H_hr.vectorL = vectorList;
if H_hr.overlap
NRPTS_S = numel(H_hr.ScoeL);
sizeScoeL = size(H_hr.ScoeL);
ScoeLtmp = reshape(H_hr.ScoeL,[NRPTS_S,1]);
[iL,jL,nL]= ind2sub(sizeScoeL,1:NRPTS_S);
vectorList_overlap = [H_hr.vectorL_overlap(nL,:),iL.',jL.'];
H_hr.ScoeL = ScoeLtmp;
H_hr.vectorL_overlap = vectorList_overlap;
end
else
H_hr.HnumL = [];
H_hr.HcoeL = sym([]);
H_hr.vectorL = [];
end
H_hr.Type = 'list';
H_hr = H_hr.simplify(options.Accuracy);
elseif  strcmp(H_hr.Type ,'list') && options.rewind
if H_hr.overlap
[vectorList_overlap,~,ic_S] = unique(H_hr.vectorL_overlap(:,1:H_hr.Dim),'rows');
NRPTS_S= size(vectorList_overlap,1);
SnumLtmp = zeros(WANNUM,WANNUM,NRPTS_S);
ScoeLtmp = sym(zeros(WANNUM,WANNUM,NRPTS_S));
end
if H_hr.num
[vectorList,~,icL] = unique(H_hr.vectorL(:,1:H_hr.Dim),'rows');
NRPTS_= size(vectorList,1);
HnumLtmp = zeros(WANNUM,WANNUM,NRPTS_);
sizemesh = [WANNUM,WANNUM,NRPTS_];
if H_hr.overlap
for n = 1:size(H_hr.vectorL_overlap,1)
SnumLtmp(H_hr.vectorL_overlap(n,H_hr.Dim+1),H_hr.vectorL_overlap(n,H_hr.Dim+2),ic_S(n)) = H_hr.SnumL(n);
end
end
else
[vectorList,~,icL] = unique(H_hr.vectorL(:,1:H_hr.Dim),'rows');
NRPTS_= size(vectorList,1);
HcoeLtmp = sym(zeros(WANNUM,WANNUM,NRPTS_));
sizemesh = [WANNUM,WANNUM,NRPTS_];
if H_hr.overlap
for n = 1:NRPTS_S
ScoeLtmp(H_hr.vectorL_overlap(n,H_hr.Dim+1),H_hr.vectorL_overlap(n,H_hr.Dim+2),ic_S(n)) = H_hr.ScoeL(n);
end
end
end
if H_hr.num
iL = double(H_hr.vectorL(:,H_hr.Dim+1));
jL = double(H_hr.vectorL(:,H_hr.Dim+2));
indL = sub2ind(sizemesh,iL,jL,icL);
HnumLtmp(indL) = H_hr.HnumL;
H_hr.HnumL = HnumLtmp;
H_hr.vectorL = vectorList;
if H_hr.overlap
H_hr.SnumL = SnumLtmp;
H_hr.vectorL_overlap = (vectorList_overlap);
end
end
if H_hr.coe
iL = double(H_hr.vectorL(:,H_hr.Dim+1));
jL = double(H_hr.vectorL(:,H_hr.Dim+2));
indL = sub2ind(sizemesh,iL,jL,icL);
HcoeLtmp(indL) = H_hr.HcoeL;
H_hr.HcoeL = HcoeLtmp;
H_hr.vectorL = vectorList;
if H_hr.overlap
H_hr.ScoeL = ScoeLtmp;
H_hr.vectorL_overlap = (vectorList_overlap);
end
end
H_hr.Type = 'mat';
else
end
end
