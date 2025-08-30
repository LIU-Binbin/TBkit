%[text] ## BlackBox
HRobj=HR.from_wannier90('wannier90_hr.dat');%里面是hr.dat的文件名 
HRobj = HRobj.input_Rm('POSCAR')
[HRobj.orbL, HRobj.elementL, HRobj.quantumL] = wout_read('wannier90.wout', 'POSCAR');
Mz = Oper.mirror([0,0,1],'TBkitObj',HRobj,'t',[0,0,0]); % check t
[HR_list,eig_list]=subsystem_construct(HRobj,Mz) %根据Mz对称性,将HR模型分成两个对应不同Mz本征值的子系统HR
HR_plus=HR_list{1}; % Mz=i的子系统HR
HR_minus=HR_list{2}; % Mz=-i的子系统HR
%[text] ### Mz子系统的BC和Chern，即总系统的MirrorBC和MirrorChern
%[text] 'BAND\_index'为占据带数目:
%[text] 例如此教程中
%[text] 总系统占据 28 带，两个 Mz 子系统各占据 14 带，所以计算 MirrorBC 使用 1:14 
%[text] #### Mz=i的子系统
[BCCAR_plus,Grid] =HR_plus.BC_2D('knum1',101,'knum2',101,'BAND_index',1:14); %根据Mz算符，计算PtCl3的Mirror-BC
BCplot2D((BCCAR_plus),Grid,double(HR_plus.Rm),'BZ',true,'BZmode','2D','ColorCut',0.05,'shading',true);
fprintf("Mirror(i) Chern number of PtCl3-AA'-AFM(y) = %f",sum(BCCAR_plus,'all')/(2*pi));
%[text] #### Mz=-i的子系统
[BCCAR_minus,Grid] =HR_minus.BC_2D('knum1',101,'knum2',101,'BAND_index',1:14); 
BCplot2D((BCCAR_minus),Grid,double(HR_minus.Rm),'BZ',true,'BZmode','2D','ColorCut',0.05,'shading',true);
fprintf("Mirror(-i) Chern number of PtCl3-AA'-AFM(y) = %f",sum(BCCAR_minus,'all')/(2*pi));
return;

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
