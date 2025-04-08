function Hcktobj = set_home(Hcktobj,Subcktobj,PortInL,PortOutL,DescriptionL)
arguments
Hcktobj Hckt;
Subcktobj Subckt;
PortInL = 1;
PortOutL = 2;
DescriptionL = Subcktobj.Description;
end
Hcktobj = Hcktobj.set_hop(zeros(1,Hcktobj.dim),Subcktobj,PortInL,PortOutL,DescriptionL);
end
