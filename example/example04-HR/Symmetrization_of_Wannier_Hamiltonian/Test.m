%% 
% Wannsymm: <https://github.com/ccao/WannSymm https://github.com/ccao/WannSymm>
% 
% TBmodels: <https://tbmodels.greschd.ch/en/latest/index.html https://tbmodels.greschd.ch/en/latest/index.html>
% 
% WannierBerri: <https://docs.wannier-berri.org/en/master/symmetries.html https://docs.wannier-berri.org/en/master/symmetries.html>
% 
% WannierTools: <http://www.wanniertools.com/symmetrization.html http://www.wanniertools.com/symmetrization.html>
%% 
% Te
% 
% 1.0
% 
% 4.4580001831         0.0000000000         0.0000000000
% 
% -2.2290000916         3.8607414086         0.0000000000
% 
% 0.0000000000         0.0000000000         5.9250001907
% 
% Te
% 
% 3
% 
% Direct
% 
% 0.745000000         0.000000000         0.333333343
% 
% 0.000000000         0.745000000         0.666666677
% 
% 0.255000000         0.255000000         0.000000000
% 
% Space Group Number: 152 Point Group: 18 [ D3 ] International: P3_121
% 
% 
% 
% DFT(vasp5.4.4 with patch: <https://github.com/Chengcheng-Xiao/VASP2WAN90_v2_fix) 
% https://github.com/Chengcheng-Xiao/VASP2WAN90_v2_fix)>
% 
% 

cd DFT/
[fig,Ax] = Figs(1,2);
EIGENCAR = EIGENVAL_read;
Ecut = [-6, 6]; 
bandplot(EIGENCAR,Ecut,'Color','r','ax',Ax(1));
%
%[PROCAR_collection,EIGENCAR] = PROCAR_read('SOC_flag',1);

[WEIGHTCAR_struct_cell, Name_list, ~] = WEIGHTCAR_read_dir('PBAND');
WEIGHTCAR_struct_Te = WEIGHTCAR_struct_cell{1}(2:4);  
titlestring = 'fatband-Te-soc'; 

% Generate band plot

[fig, ax] = pbandplot(WEIGHTCAR_struct_Te, EIGENCAR,...
    'Ecut', Ecut,...
    'title', titlestring,...
    'silent', true,'ax',Ax(2));
set(fig,'Position',[0 0 2400 1000]);
set(Ax(1),'FontSize',32);
set(Ax(2),'FontSize',32);
% Generate .win & wt.in
init_v2w(1,[-6,6],'check',1,'Projector_list',[1,1,1],'BandIndex',[7,7+6*3+1]);
selectbands(EIGENCAR,eye(3),'dis_check','BandIndex',[7,7+6*3-1])
cd ..
%% 
% wannier
% 
% 
% 
% shift FermiEnergy！/

cd Wan/
Hr_Te = HR.from_wannier90();
Hr_Te.Rm = POSCAR_read;
[Hr_Te.orbL,Hr_Te.elementL, Hr_Te.quantumL] = wout_read;
EIGENCAR_wan = Hr_Te.EIGENCAR_gen()-0.02;
[fig,Ax] = Figs(1,2);
bandplot({EIGENCAR,EIGENCAR_wan},Ecut,'Color',[0 0 1;1 0 0],'legends',["DFT","WAN"],'ax',Ax(1));
bandplot({EIGENCAR,EIGENCAR_wan},[-1,1],'Color',[0 0 1;1 0 0],'legends',["DFT","WAN"],'ax',Ax(2));
set(fig,'Position',[0 0 2400 1000]);
set(Ax(1),'FontSize',32);
set(Ax(2),'FontSize',32);


cd ..
%%
[H_soc_full, lambda_syms] = soc_term_udud_add(Hr_Te.elementL, Hr_Te.quantumL,'mode','basis');

H_soc_full
%%
BasisFunction = BasisFunc(Hr_Te)
%% 
% 

C3_  = Oper.rotation(1/3,[0,0,1],false,'t',[0 0 1/3],'sym',0)
C2_y = Oper.rotation(1/2,[1,0,0],false,'t',[0 0 2/3],'sym',0)
T = Oper.time_reversal(3)
Gen = [C3_,C2_y,T ];
G = Gen.generate_group
for i = 1:numel(Gen)
    if isnan(Gen(i).U)
        Gen(i) = Gen(i).attachRm(Hr_Te.Rm);
        try
            Gen(i).U = BasisFunction.rotation('Oper',Gen(i),'Rm',Hr_Te.Rm);
        catch
            error('fail to generate Oper matrix!');
        end
    end
end
Gen(1)
Gen(2)
det(Gen(2).U)
Gen(3)
%%
Hr_Te_list = Hr_Te.rewrite()
Hr_Te_list = Hr_Te_list.simplify
%%
cd Wan/
%%
Hr_Te_sym = Hr_Te_list.applyOper(Gen,"generator",1)
%%
EIGENCAR_wan2 = Hr_Te_sym.EIGENCAR_gen()-0.02;
[fig,Ax] = Figs(1,2);
bandplot({EIGENCAR,EIGENCAR_wan},[-1,1],'Color',[0 0 1;1 0 0],'legends',["DFT","WAN"],'ax',Ax(1),'title','Before');
bandplot({EIGENCAR,EIGENCAR_wan2},[-1,1],'Color',[0 0 1;1 0 0],'legends',["DFT","WAN"],'ax',Ax(2),'title','After: Sym By TBkit');
set(fig,'Position',[0 0 2400 1000]);
set(Ax(1),'FontSize',32);
set(Ax(2),'FontSize',32);