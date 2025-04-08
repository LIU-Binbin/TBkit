function H_hr = addorb(H_hr,orblist,options)
arguments
H_hr HR;
orblist;
options.inner  = true;
end
WANNUM = H_hr.WAN_NUM;
norb = length(orblist);
HnumLtmp = H_hr.HnumL;
HnumLtmp2 = zeros(size(H_hr.HnumL)+[norb,norb,0]);
HnumLtmp2(1:WANNUM,1:WANNUM,:) = HnumLtmp;
HnumLtmp2(WANNUM+1:WANNUM+norb,1:WANNUM,:) = HnumLtmp(orblist,:,:);
HnumLtmp2(1:WANNUM,WANNUM+1:WANNUM+norb,:) = HnumLtmp(:,orblist,:);
H_hr.HnumL = HnumLtmp2;
end
