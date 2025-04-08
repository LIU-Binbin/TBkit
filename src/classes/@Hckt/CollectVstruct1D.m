function [DOSCAR_2D,klist1,OmgL] = CollectVstruct1D(ObservationsMat,VectorList,OmegaCut,SpectrumX)
if nargin < 2
OmegaCar = ObservationsMat.';
EIGENCAR = abs(fft((OmegaCar),[],2));
DOSCAR_2D = EIGENCAR(:,[1:end,1]);
klist1 = [];
OmgL = [];
else
VectorList1 =  VectorList;
ObservationsMat = VectorList;
VectorList = ObservationsMat;
OmegaCar = ObservationsMat.';
EIGENCAR = abs(fft((OmegaCar),[],2));
DOSCAR_2D = EIGENCAR(:,[1:end,1]);
klist1 = [];
OmgL = [];
end
end
