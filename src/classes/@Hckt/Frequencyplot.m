function [fig,ax] = Frequencyplot(ObservationsMat,SpectrumL,OmegaCut,options)
arguments
ObservationsMat;
SpectrumL;
OmegaCut = [0,2e4];
options.fig =  handle([]);
options.ax =  handle([]);
options.FontName = 'Helvetica';
options.FontSize = 24;
options.Units = 'pixels';
options.Position = [];
options.Color =  [rand rand rand];
options.title = '';
options.xlabel='Frequency (Hz)';
options.ylabel='';
end
if isempty(options.fig) && isempty(options.ax)
[fig,ax] = TBkit_tool.creat_figure('Units',options.Units,'Position',options.Position);
else
fig = options.fig;
ax = options.ax;
end
for i = 1:size(ObservationsMat,1)
plot(SpectrumL,ObservationsMat(i,:),'DisplayName',num2str(i),'Color',options.Color);
hold on
end
title(ax,options.title);
xlabel(ax,options.xlabel);
ylabel(ax,options.ylabel);
xlim(OmegaCut);
end
