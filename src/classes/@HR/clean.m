function H_hr = clean(H_hr,WANNUM)
if nargin < 2
WANNUM = H_hr.WAN_NUM;
end
if strcmp(H_hr.Type,'Sparse')
elseif strcmp(H_hr.Type,'mat')
H_hr.HnumL = zeros(WANNUM,WANNUM,H_hr.NRPTS);
if H_hr.coe
H_hr.HcoeL = sym(H_hr.HnumL);
else
H_hr.HcoeL =sym([]);
end
if H_hr.overlap
H_hr.SnumL = zeros(WANNUM,WANNUM,H_hr.NRPTS);
if H_hr.coe
H_hr.ScoeL = sym(H_hr.SnumL);
else
H_hr.ScoeL =sym([]);
end
end
end
end
