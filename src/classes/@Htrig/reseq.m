function H_htrig = reseq(H_htrig,wan_list,kinds_list)
if nargin < 3
kinds_list = ':';
end
if ~isequal(wan_list,':')
switch H_htrig.Type
case 'sparse'
for i = 1:H_htrig.Kinds
H_htrig.HnumL{i}=H_htrig.HnumL{i}(wan_list,wan_list);
end
case 'list'
wan_list =  all(ismember(H_htrig.vectorL(:,H_htrig.Dim+1:H_htrig.Dim+2),int32(wan_list)),2);
if H_htrig.num
H_htrig.HnumL=H_htrig.HnumL(wan_list,:);
H_htrig.HsymL_numL=H_htrig.H_htrig.HsymL_numL(wan_list,:);
H_htrig.HsymL_coeL = [];
H_htrig.HcoeL= [];
end
if H_htrig.coe
H_htrig.HcoeL=H_htrig.HcoeL(wan_list,:);
H_htrig.HsymL_coeL=H_htrig.H_htrig.HsymL_coeL(wan_list,:);
H_htrig.HsymL_numL = [];
H_htrig.HnumL= [];
end
case 'mat'
if H_htrig.num
H_htrig.HnumL=H_htrig.HnumL(wan_list,wan_list,:);
else
H_htrig.HnumL = [];
H_htrig.HsymL_numL = [];
end
if H_htrig.coe
H_htrig.HcoeL=H_htrig.HcoeL(wan_list,wan_list,:);
else
H_htrig.HcoeL= [];
H_htrig.HcoeL = sym([]);
end
otherwise
if H_htrig.num
H_htrig.HnumL=H_htrig.HnumL(wan_list,wan_list,:);
else
H_htrig.HnumL = [];
end
if H_htrig.coe
H_htrig.HcoeL=H_htrig.HcoeL(wan_list,wan_list,:);
else
H_htrig.HcoeL = sym([]);
end
end
if ~isempty(H_htrig.sites)
try
H_htrig.sites = H_htrig.sites(wan_list);
catch
end
end
if ~isempty( H_htrig.orbL )
H_htrig.orbL = H_htrig.orbL(wan_list,:);
end
if ~isempty( H_htrig.elementL )
H_htrig.elementL = H_htrig.elementL(wan_list,:);
end
if ~isempty( H_htrig.quantumL)
H_htrig.quantumL = H_htrig.quantumL(wan_list,:);
end
end
if ~isequal(kinds_list,':')
switch H_htrig.Type
case 'sparse'
H_htrig.HnumL=H_htrig.HnumL(kinds_list);
case 'list'
if H_htrig.num
H_htrig.HsymL_numL=H_htrig.HsymL_numL(kinds_list,:);
H_htrig.HnumL=H_htrig.HnumL(kinds_list);
if ~isempty(H_htrig.HcoeL)
try
H_htrig.HcoeL=H_htrig.HcoeL(:,:,kinds_list);
catch
end
end
if ~isempty(H_htrig.HsymL_coeL)
try
H_htrig.HsymL_coeL=H_htrig.HsymL_coeL(kinds_list,:);
catch
end
end
end
if H_htrig.coe
H_htrig.HsymL_coeL=H_htrig.HsymL_coeL(kinds_list,:);
H_htrig.HcoeL=H_htrig.HcoeL(kinds_list);
if ~isempty(H_htrig.HnumL)
try
H_htrig.HnumL=H_htrig.HnumL(:,:,kinds_list);
catch
end
end
if ~isempty(H_htrig.HsymL_numL)
try
H_htrig.HsymL_numL=H_htrig.HsymL_numL(kinds_list,:);
catch
end
end
end
case  'mat'
if H_htrig.num
H_htrig.HsymL_numL=H_htrig.HsymL_numL(kinds_list,:);
H_htrig.HnumL=H_htrig.HnumL(:,:,kinds_list);
if ~isempty(H_htrig.HcoeL)
try
H_htrig.HcoeL=H_htrig.HcoeL(:,:,kinds_list);
catch
end
end
if ~isempty(H_htrig.HsymL_coeL)
try
H_htrig.HsymL_coeL=H_htrig.HsymL_coeL(kinds_list,:);
catch
end
end
end
if H_htrig.coe
H_htrig.HsymL_coeL=H_htrig.HsymL_coeL(kinds_list,:);
H_htrig.HcoeL=H_htrig.HcoeL(:,:,kinds_list);
if ~isempty(H_htrig.HnumL)
try
H_htrig.HnumL=H_htrig.HnumL(:,:,kinds_list);
catch
end
end
if ~isempty(H_htrig.HsymL_numL)
try
H_htrig.HsymL_numL=H_htrig.HsymL_numL(kinds_list,:);
catch
end
end
end
otherwise
H_htrig.HsymL_trig = H_htrig.HsymL_trig(kinds_list);
if H_htrig.coe
H_htrig.HcoeL = H_htrig.HcoeL(:,:,kinds_list);
end
if H_htrig.num
H_htrig.HnumL = H_htrig.HnumL(:,:,kinds_list);
end
end
end
end
