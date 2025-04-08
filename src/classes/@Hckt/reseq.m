function HcktObj = reseq(HcktObj,seqL)
HcktObj.ScktL = HcktObj.ScktL(seqL,:);
HcktObj.vectorL = HcktObj.vectorL(seqL);
HcktObj.PortInCell = HcktObj.PortInCell(seqL,:);
HcktObj.PortOutCell = HcktObj.PortOutCell(seqL,:);
HcktObj.DescriptionL = HcktObj.DescriptionL(seqL,:);
end
