function H_hr = from_Hsparse(Hsparse)
vectorL = Hsparse.vectorL;
[NRPTS,~] = size(vectorL);
WAN_NUM = length(Hsparse.HnumL{1});
HcoeL = [];
HnumL = zeros(WAN_NUM,WAN_NUM,NRPTS);
for i = 1:NRPTS
HnumL(:,:,i) = full(Hsparse.HnumL{i});
end
H_hr = HR(WAN_NUM,vectorL,'HnumL',HnumL,'HcoeL',HcoeL);
end
