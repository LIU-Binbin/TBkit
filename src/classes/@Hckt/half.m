function HcktObj = half(HcktObj,checkcommute)
if nargin < 2
checkcommute = false;
end
MaxR = max(abs(HcktObj.vectorAll),[],'all');
vectorLcheck = Hckt.dim2vectorL(HcktObj.dim,MaxR);
[whichexist,~] = ismember(HcktObj.vectorAll ,vectorLcheck,'rows') ;
ToKeepRef = find(whichexist);
ToKeepL = ismember(HcktObj.vectorL,ToKeepRef);
if checkcommute
commuteL = [HcktObj.ScktL.commute];
ToKeepL = logical(ToKeepL + ~commuteL);
end
HcktObj = HcktObj.reseq(ToKeepL);
end
