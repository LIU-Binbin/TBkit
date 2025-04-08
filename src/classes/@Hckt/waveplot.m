function [ax,WaveFunc] = waveplot(ObservationsMat,VectorList,SpectrumL,orbL,options_select,options)
arguments
ObservationsMat = [];
VectorList = [];
SpectrumL = [];
orbL = [];
options_select.Frequency = -1;
options_select.Width = 0;
options_select.scale = 1;
options.ax =  handle([]);
options.Rm = [];
options.POSCAR = 'POSCAR';
options.WaveMin = 1e-3;
options.WaveColor = 'r';
options.WaveSize = 1;
options.OrbColor = 'k';
options.OrbSize = 1;
end
optionscell = namedargs2cell(options);
if options_select.Frequency == -1
options_select.Frequency = SpectrumL(end/2);
end
if options_select.Width == 0
Scale = round(log(options_select.Frequency)/log(10));
options_select.Width = 10^(Scale-3);
end
OmegaCut = [options_select.Frequency - options_select.Width,options_select.Frequency + options_select.Width];
ChooseL = SpectrumL >= OmegaCut(1) & SpectrumL <= OmegaCut(2);
[VectorList,ObservationsMat] = HollowKnight.generalcontractrow2(VectorList,ObservationsMat);
ObservationsMat = abs(ObservationsMat);
NVectorList = size(VectorList,1);
RvectorL = VectorList(:,1:end-1);
NodeL = VectorList(:,end);
Dim = size(RvectorL,2);
if Dim <3
RvectorL = [RvectorL,zeros(NVectorList,3-Dim)];
end
littleorbL = orbL(NodeL,:);
ORBL = RvectorL+littleorbL;
WaveFunc = normalize(sum(ObservationsMat(:,ChooseL),2),'range',[0,1])*options_select.scale;
ax = waveplot(ORBL,WaveFunc,optionscell{:});
end
