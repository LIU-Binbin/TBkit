function [fig,ax] = show(H_hr,mode,options)
arguments
H_hr HR;
mode char {mustBeMember(mode,{'POSCAR','NRPTS','HOPPING'})} = 'HOPPING';
options.fig = figure('WindowState','maximized');
options.scale = 1;
options.atomscale = 0.3;
options.TwoD = false;
options.vectorList = [];
options.fast = true;
options.ax = [];
options.Select = [];
options.cmap = @parula;
options.Title = '';
end
Rm_ = H_hr.Rm*options.scale;
if H_hr.vectorhopping
H_hr = H_hr.GenfromOrth();
end
if H_hr.Dim > 3
error('Dont support Dim >3');
end
H_hr.orbL = double(H_hr.orbL);
switch mode
case 'POSCAR'
fig  = POSCAR_plot(Rm_,H_hr.sites,H_hr.Atom_name,H_hr.Atom_num);
case 'NRPTS'
if isempty(options.vectorList)
vectorList = double(unique(H_hr.vectorL(:,1:H_hr.Dim),'rows'));
else
vectorList = options.vectorList ;
end
NRPTS_ = size(vectorList,1);
positions = H_hr.orbL;
sites_ = [];
if isempty(H_hr.elementL)
elementList = 6*ones(size(positions,1),1);
else
elementList = H_hr.elementL;
end
Atom_name_ = repmat(elementList,[NRPTS_,1]);
for i = 1:NRPTS_
sites_ = [sites_;positions+vectorList(i,:)] ;
end
[fig,ax]  = POSCAR_plot(Rm_,sites_,Atom_name_,...
'vectorL',vectorList,...
'fig',options.fig,...
'ax',options.ax,...
'scale',options.scale,...
'atomscale',options.atomscale,...
'TwoD',options.TwoD ,...
'fast',options.fast);
case 'HOPPING'
H_hr = H_hr.rewrite();
[fig,ax] = H_hr.show('NRPTS',...
'scale',options.scale,...
'ax',options.ax,...
'atomscale',options.atomscale,...
'fig',options.fig,...
'TwoD',options.TwoD );
Rm_ = double(Rm_ * options.scale);
if H_hr.coe
if isempty(options.Select)
Plot_Hopping(H_hr.vectorL,H_hr.HcoeL,Rm_,H_hr.orbL,'ax',ax,'TwoD',options.TwoD,'cmap', options.cmap );
else
SelectL = ismember(H_hr.HcoeL,options.Select);
Plot_Hopping(H_hr.vectorL(SelectL,:),H_hr.HcoeL(SelectL),Rm_,H_hr.orbL,'ax',ax,'TwoD',options.TwoD ,'cmap', options.cmap);
end
else
if isempty(options.Select)
Plot_Hopping(H_hr.vectorL,H_hr.HnumL,Rm_,H_hr.orbL,'ax',ax,'TwoD',options.TwoD );
else
SelectL = H_hr.HnumL == options.Select;
Plot_Hopping(H_hr.vectorL(SelectL,:),H_hr.HnumL(SelectL),Rm_,H_hr.orbL,'ax',ax,'TwoD',options.TwoD );
end
end
if options.TwoD
view(ax,0,90);
l = light(ax);
l.Position = [0,0,1];
end
otherwise
end
title(options.ax,options.Title);
end
