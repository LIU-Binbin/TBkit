%Spin Polarized bandplot for Altermagnet Ca(CoN)2
%you need prepare POSCAR,KPOINTS & EIGENVAL from band calculatrion,and DOSCAR from scf 
Efermi=GetFermi('vasp','DOSCAR_CaCo2N2');  %get Fermi energy from DOSCAR
%When use spin polarized(ISPIN=2),you got 2 columns of EIGENVAL
[EIGENCAR1,EIGENCAR] = EIGENVAL_read("vasp",'EIGENVAL_CaCo2N2',Efermi); %read energy eigen value
bandplot({EIGENCAR1,EIGENCAR},[-1,1],'KPOINTS','KPOINTS_CaCo2N2','POSCAR','POSCAR_CaCo2N2',...
    'Color',[0 0 1;1 0 0],'LineSpec', {'--','-'})  
title('Ca(CoN)2-SpinPolarized')