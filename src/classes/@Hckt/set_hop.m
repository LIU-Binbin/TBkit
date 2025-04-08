function Hcktobj = set_hop(Hcktobj,vector,Subcktobj,PortInL,PortOutL,DescriptionL,options)
arguments
Hcktobj Hckt;
vector;
Subcktobj Subckt;
PortInL = 1;
PortOutL = 2;
DescriptionL = Subcktobj.Description;
options.Port2PortMode = false;
end
[~,vectorlabel] = ismember(vector,Hcktobj.vectorAll,'rows');
if vectorlabel == 0
Hcktobj = add_empty_one(Hcktobj,vector);
vectorlabel = size(Hcktobj.vectorAll,1);
end
if vectorlabel == 1
seq =1;
else
seq = Hcktobj.NRPTS+1;
end
Hcktobj.ScktL(seq,1) = Subcktobj;
Hcktobj.vectorL(seq) = vectorlabel;
Hcktobj.PortInCell{seq,1} = PortInL;
Hcktobj.PortOutCell{seq,1}= PortOutL;
Hcktobj.DescriptionL{seq,1} = char(DescriptionL);
if options.Port2PortMode
checkvectorL = Hckt.dim2vectorL(Hcktobj.dim);
if ~all(vector==0)
if ismember(vector,checkvectorL,'rows')
Hcktobj.ScktL(seq,1).ReverseConnection = false;
Hcktobj.Port2VectorDist(PortInL) = vector;
Hcktobj.Port2PortDist(PortInL) = PortOutL;
Hcktobj.Port2PortDistForVectorise(PortInL) = PortOutL;
Hcktobj.Port2PortDist_BinBin(PortInL) = PortInL;
else
Hcktobj.ScktL(seq,1).ReverseConnection = true;
Hcktobj.Port2VectorDist(PortInL) =  zeros(1,Hcktobj.dim);
Hcktobj.Port2PortDist(PortInL) = PortInL;
Hcktobj.Port2PortDistForVectorise(PortInL) = PortOutL;
Hcktobj.Port2PortDist_BinBin(PortInL) = PortOutL;
end
end
end
end
