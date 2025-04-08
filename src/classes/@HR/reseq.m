function H_hr = reseq(H_hr,wan_list,nrpt_list,nrpt_list_S)
if nargin < 3
nrpt_list = ':';
end
if nargin < 4
nrpt_list_S = ':';
end
if ~isequal(wan_list,':')
if strcmp(H_hr.Type,'sparse')
for i = 1:H_hr.NRPTS
H_hr.HnumL{i}=H_hr.HnumL{i}(wan_list,wan_list);
end
elseif strcmp(H_hr.Type,'list')
wan_list =  all(ismember(H_hr.vectorL(:,H_hr.Dim+1:H_hr.Dim+2),(wan_list)),2);
if H_hr.num
H_hr.HnumL=H_hr.HnumL(wan_list,:);
end
if H_hr.coe
H_hr.HcoeL=H_hr.HcoeL(wan_list,:);
end
H_hr.vectorL=H_hr.vectorL(wan_list,:);
if H_hr.overlap
H_hr.SnumL=H_hr.SnumL(wan_list,:);
H_hr.ScoeL=H_hr.ScoeL(wan_list,:);
H_hr.vectorL_overlap=H_hr.vectorL_overlap(wan_list,:);
end
elseif strcmp(H_hr.Type,'mat')
if H_hr.num
H_hr.HnumL=H_hr.HnumL(wan_list,wan_list,:);
if H_hr.overlap
H_hr.SnumL=H_hr.SnumL(wan_list,wan_list,:);
end
else
H_hr.HnumL = [];
if H_hr.overlap
H_hr.SnumL=[];
end
end
if H_hr.coe
H_hr.HcoeL=H_hr.HcoeL(wan_list,wan_list,:);
if H_hr.overlap
H_hr.ScoeL=H_hr.ScoeL(wan_list,wan_list,:);
end
else
H_hr.HcoeL = sym([]);
if H_hr.overlap
H_hr.ScoeL=sym([]);
end
end
end
if ~isempty(H_hr.sites)
try
H_hr.sites = H_hr.sites(wan_list);
catch
end
end
if ~isempty( H_hr.orbL )
H_hr.orbL = H_hr.orbL(wan_list,:);
end
if ~isempty( H_hr.elementL )
H_hr.elementL = H_hr.elementL(wan_list,:);
end
if ~isempty( H_hr.quantumL)
H_hr.quantumL = H_hr.quantumL(wan_list,:);
end
end
if ~isequal(nrpt_list,':')
H_hr.vectorL=H_hr.vectorL(nrpt_list,:);
if H_hr.overlap
H_hr.vectorL_overlap=H_hr.vectorL_overlap(nrpt_list_S,:);
end
if strcmp(H_hr.Type,'sparse')
H_hr.HnumL=H_hr.HnumL(nrpt_list);
elseif strcmp(H_hr.Type,'list')
if H_hr.vectorhopping
H_hr.AvectorL=H_hr.AvectorL(nrpt_list,:);
H_hr.BvectorL=H_hr.BvectorL(nrpt_list,:);
CL1 = H_hr.CvectorL(1:end/2,:);
CL2 = H_hr.CvectorL(end/2+1:end,:);
H_hr.CvectorL = [CL1(nrpt_list,:);CL2(nrpt_list,:)];
return;
end
if H_hr.num
H_hr.HnumL=H_hr.HnumL(nrpt_list);
end
if H_hr.coe
H_hr.HcoeL=H_hr.HcoeL(nrpt_list);
end
if H_hr.overlap
if H_hr.num
H_hr.SnumL=H_hr.SnumL(nrpt_list_S,:);
end
if H_hr.coe
H_hr.ScoeL=H_hr.ScoeL(nrpt_list_S,:);
end
end
elseif strcmp(H_hr.Type,'mat')
if H_hr.num
H_hr.HnumL=H_hr.HnumL(:,:,nrpt_list);
if H_hr.overlap
H_hr.SnumL=H_hr.SnumL(:,:,nrpt_list);
end
end
if H_hr.coe
H_hr.HcoeL=H_hr.HcoeL(:,:,nrpt_list);
if H_hr.overlap
H_hr.ScoeL=H_hr.ScoeL(:,:,nrpt_list);
end
end
end
end
end
