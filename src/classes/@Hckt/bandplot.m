function [ax] = bandplot(EIGNECAR,OmegaCut,SpectrumL,options,optionskpath)
arguments
EIGNECAR = [];
OmegaCut = [0 20000];
SpectrumL = [];
options.ax =  handle([]);
options.FontName = 'Helvetica';
options.FontSize = 24;
options.Units = 'pixels';
options.Position = [];
options.Color =  parula;
options.title = '';
options.xlabel='';
options.ylabel='Frequency (Hz)';
options.alpha = 1;
options.KPOINTS = 'KPOINTS';
options.POSCAR = 'POSCAR';
options.shift = false;
optionskpath.klist_l = [];
optionskpath.kpoints_l = [];
optionskpath.kpoints_name = [];
end
if isempty(options.ax)
Fig = create_figure('Position',options.Position);
ax = Fig.axes;
else
ax = options.ax;
end
if exist(options.KPOINTS,'file') && exist(options.POSCAR,'file') && isempty(optionskpath.kpoints_name )
Rm=POSCAR_read(options.POSCAR);
[~,klist_l,~,kpoints_l,kpoints_name]=kpathgen3D(Rm,options.KPOINTS);
UseKpath = true;
if options.shift
X = linspace(min(klist_l),max(klist_l),size(EIGNECAR,2)+1);
else
X = klist_l;
end
elseif ~isempty(optionskpath.kpoints_name) && ~isempty(optionskpath.kpoints_l)
klist_l = optionskpath.klist_l;
kpoints_name = optionskpath.kpoints_name;
kpoints_l = optionskpath.kpoints_l;
UseKpath = true;
if options.shift
X = linspace(min(klist_l),max(klist_l),size(EIGNECAR,2)+1);
else
X = klist_l;
end
else
UseKpath = false;
if options.shift
X = [1:size(EIGNECAR,2),size(EIGNECAR,2)+1];
else
X = [1:size(EIGNECAR,2),size(EIGNECAR,2)];
end
end
ChooseL = SpectrumL>=OmegaCut(1) & SpectrumL<=OmegaCut(2);
Y = SpectrumL(ChooseL);
if options.shift
Z = EIGNECAR(ChooseL,[1:size(EIGNECAR,2),1]);
else
Z = EIGNECAR(ChooseL,:);
end
[Xgrid,Ygrid] = meshgrid(X,Y);
surf(ax,Xgrid,Ygrid,Z,'EdgeColor','none','FaceAlpha',options.alpha);
axis(ax,'tight');
view(ax,2);
colormap(options.Color);
title(ax,options.title);
if UseKpath
Xcut = [min(X),max(X)];
Ycut = [min(Y),max(Y)];
ax= set_reference(kpoints_l,kpoints_name,Xcut,Ycut,...
'ax',ax,...
'xlabel',options.xlabel,...
'ylabel',options.ylabel ...
);
else
xlabel(ax,options.xlabel);
ylabel(ax,options.ylabel);
end
end
