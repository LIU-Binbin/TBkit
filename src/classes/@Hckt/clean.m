function HcktObj = clean(HcktObj)
hollowlabel = find(HcktObj.ScktL == 0);
HcktObj = HcktObj.reseq(hollowlabel);
end
