function H_hk  = init(H_hk)
sizeH = size(H_hk.HcoeL);
HcoeL_tmp = sym('A',sizeH,'real')+1i*sym('B',sizeH,'real');
for i =1:H_hk.Kinds
HcoeL_tmp(:,:,i) = triu(HcoeL_tmp(:,:,i))+triu(HcoeL_tmp(:,:,i),1)';
end
H_hk.HcoeL = HcoeL_tmp;
end
