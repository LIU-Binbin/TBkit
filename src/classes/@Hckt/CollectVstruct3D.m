function [DOSCAR_4D,klist1,klist2,klist3,OmgL] = CollectVstruct3D(VectorList,ObservationsMat,OmegaCut,SpectrumX,options)
arguments
VectorList
ObservationsMat
OmegaCut
SpectrumX
options.fin_dir = [0,0,0];
end
SelectOmegaL = find(SpectrumX >= OmegaCut(1) & SpectrumX <= OmegaCut(2));
Rvector = VectorList(:,1:end-1);
if options.fin_dir(1)
NX = 1;
klist1 = 1;
else
NX = max(Rvector(:,1));
klist1 = 1:NX+1;
end
if options.fin_dir(2)
NY = 1;
klist2 = 1
else
NY = max(Rvector(:,2));
klist2 = 1:NY+1;
end
if options.fin_dir(2)
NZ = 1;
klist3 = 1;
else
NZ = max(Rvector(:,3));
klist3 = 1:NZ+1;
end
sizemesh = [NX,NY,NZ];
OmgL = SpectrumX(SelectOmegaL);
ObservationsMat = ObservationsMat(:,SelectOmegaL);
DOSCAR_4D = zeros(sizemesh(1),sizemesh(2),sizemesh(3),length(OmgL));
for i = 1:numel(SelectOmegaL)
tmpDATA = abs(fftn(reshape(ObservationsMat(:,i),NZ,NY,NX)));
DOSCAR_4D(:,:,:,i) = permute(tmpDATA,[3,2,1]);
end
if ~options.fin_dir(1)
DOSCAR_4D(sizemesh(1)+1,:,:,:) = DOSCAR_4D(1,:,:,:);
end
if ~options.fin_dir(2)
DOSCAR_4D(:,sizemesh(2)+1,:,:) = DOSCAR_4D(:,1,:,:);
end
if ~options.fin_dir(3)
DOSCAR_4D(:,:,sizemesh(3)+1,:) = DOSCAR_4D(:,:,1,:);
end
end
