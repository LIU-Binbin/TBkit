function [DOSCAR_3D,klist1,klist2,OmgL] = CollectVstruct2D(VectorList,ObservationsMat,OmegaCut,SpectrumX)
SelectOmegaL = find(SpectrumX >= OmegaCut(1) & SpectrumX <= OmegaCut(2));
Rvector = VectorList(:,1:end-1);
NX = max(Rvector(:,1));
NY = max(Rvector(:,2));
sizemesh = [NX,NY];
klist1 = 1:NX+1;
klist2 = 1:NY+1;
OmgL = SpectrumX(SelectOmegaL);
ObservationsMat = ObservationsMat(:,SelectOmegaL);
DOSCAR_3D = zeros(sizemesh(1),sizemesh(2),length(OmgL));
for i = 1:numel(SelectOmegaL)
DOSCAR_3D(:,:,i) = abs(fft2(reshape(ObservationsMat(:,i),NX,NY)));
end
DOSCAR_3D(sizemesh(1)+1,:,:) = DOSCAR_3D(1,:,:);
DOSCAR_3D(:,sizemesh(2)+1,:) = DOSCAR_3D(:,1,:);
end
