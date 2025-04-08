function H_hr = from_Hstruct(Hstruct)
[NRPTS,~]=size(Hstruct);
try
WAN_NUM = length(Hstruct(1).Hnum);
catch
WAN_NUM = length(Hstruct(1).Hcoe);
end
V = [Hstruct.vector];
vectorL =(reshape(V,3,length(V)/3)');
HnumL = zeros(WAN_NUM,WAN_NUM,NRPTS);
HcoeL = sym(zeros(WAN_NUM,WAN_NUM,NRPTS));
for i = 1:NRPTS
try
HnumL(:,:,i) = Hstruct(i).Hnum/Hstruct(i).Degen;
catch
HnumL(:,:,i) = zeros(WAN_NUM)  ;
end
try
HcoeL(:,:,i) = Hstruct(i).Hcoe;
catch
end
end
H_hr = HR(WAN_NUM,vectorL,'HnumL',HnumL,'HcoeL',HcoeL);
end
