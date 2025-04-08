function [VectorList,ObservationsMat,SpectrumX] = FakeSimulation(H_hr,meshs,OmegaCut,options)
arguments
H_hr HR;
meshs double;
OmegaCut double;
options.CktDim = 3;
options.Dim = 3 ;
options.wnum = 1000 ;
options.eta = 1e-6  ;
options.EIGENFACTOR =1E-15;
options.fin_dir = [];
end
if isempty(options.fin_dir)
options.fin_dir = zeros(1,options.Dim);
end
wnum = options.wnum;
if OmegaCut(1) == 0
OmegaCut(1) = 0.1*10^round(log(OmegaCut(2))/log(10));
end
Ecut = sort((OmegaCut*2*pi).^(-2)/(options.EIGENFACTOR));
SpectrumX = linspace(OmegaCut(1),OmegaCut(2),options.wnum);
SpectrumE = linspace(Ecut(1),Ecut(2),options.wnum);
H_hr = H_hr.ForceToMat();
if options.CktDim ~= length(meshs)
options.CktDim = length(meshs);
end
if  options.Dim ~= H_hr.Dim
options.Dim = H_hr.Dim;
end
if options.CktDim ~= options.Dim
meshs = [meshs,ones(1,options.Dim-options.CktDim)];
end
Ns = diag(meshs);
H_hr_super = H_hr.supercell_hr(Ns,"OBC",options.fin_dir);
sc_orb = H_hr_super.orbL;
VectorList = floor(sc_orb.*meshs);
VectorList = [VectorList(:,1:options.CktDim)+1,repmat((1:H_hr.WAN_NUM).',[prod(meshs) 1])];
nV = size(VectorList,1);
ObservationsMat = zeros(nV,wnum);
TargetH = sum(H_hr_super.HnumL(:,:,3));
Ieye = eye(nV);
for i = 1:wnum
ObservationsMat(:,i) = -imag(diag(Ieye/(Ieye*(SpectrumE(i)+1i*options.eta)-TargetH)));
end
end
