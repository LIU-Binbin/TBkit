function [DOSCAR_ND,klistCell,OmgL] = CollectVstruct(VectorList,ObservationsMat,OmegaCut,SpectrumX,options)
arguments
VectorList
ObservationsMat
OmegaCut
SpectrumX
options.fin_dir = [0,0,0];
options.cellmode = false;
end
SelectOmegaL = find(SpectrumX >= OmegaCut(1) & SpectrumX <= OmegaCut(2));
Rvector = VectorList(:,1:end-1);
Dim = size(Rvector,2);
if Dim ~= size(options.fin_dir,2)
findir = zeros(1,Dim);
else
findir = options.fin_dir;
end
if Dim < 4 && ~options.cellmode
switch Dim
case 1
[DOSCAR_ND,klist1,OmgL] = Hckt.CollectVstruct1D(ObservationsMat,VectorList,OmegaCut,SpectrumX);
klistCell{1} = klist1;
case 2
[DOSCAR_ND,klist1,klist2,OmgL] = Hckt.CollectVstruct2D(VectorList,ObservationsMat,OmegaCut,SpectrumX);
klistCell{1} = klist1; klistCell{2} = klist2;
case 3
[DOSCAR_ND,klist1,klist2,klist3,OmgL] = Hckt.CollectVstruct3D(VectorList,ObservationsMat,OmegaCut,SpectrumX);
klistCell{1} = klist1; klistCell{2} = klist2; klistCell{3} = klist3;
end
return;
end
count = Dim;
ND{count}=1;
klistCell{count} = 1;
for LNxyzw = flip(logical(findir))
if LNxyzw
ND{count} = 1;
klistCell{count} = 1;
else
ND{count}  = max(Rvector(:,count));
klistCell{count} = 1:ND{count}+1;
end
count = count - 1;
end
sizemesh = fold(@horzcat,ND);
flipND       = flip(ND);
permuteList = Dim:-1:1;
OmgL = SpectrumX(SelectOmegaL);
ObservationsMat = ObservationsMat(:,SelectOmegaL);
DOSCAR_ND{length(OmgL)} = zeros(ND{:});
for i = 1:numel(SelectOmegaL)
tmpDATA = abs(fftn(reshape(ObservationsMat(:,i),flipND{:})));
DOSCAR_ND{i} = permute(tmpDATA,permuteList);
end
end
