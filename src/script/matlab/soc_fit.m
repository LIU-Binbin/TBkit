%% LNONCOLLINEAR = T, num_iter = 1000
Ef_nsoc = 5.6944;
hr_nsoc = HR.from_wannier90('wannier90_hr.dat');
hr_nsoc = hr_nsoc.shift_Fermi_energy(Ef_nsoc);
%% LSORBIT = T
Ef_soc = 5.6974;
EIG_DFT_soc = EIGENVAL_read("vasp",'EIGENVAL_soc');
EIG_DFT_soc = EIG_DFT_soc - Ef_soc;
%% plot
EIG_hr = hr_nsoc.EIGENCAR_gen();

bandplot({EIG_DFT_soc, EIG_hr},'Color', [1 0 0; 0 0 1]);
legend(["DFT", "wannier"])
title("before fitting")
%% on-site SOC term
[H_soc_sym_udud, lambda_syms] = hr_nsoc.SOC_on_site_gen;
H_soc_sym = H_soc_sym_udud;
%% band and kpoints range
upper_bound = 24;
lower_bound = 1;
exclude_bands = 0;
options_extra.NBAND_range_DFT =(lower_bound:upper_bound) + exclude_bands;
options_extra.NBAND_range     = lower_bound:upper_bound;
options_extra.E_range = [];

options_extra.klist_range = ':';
%% search
% 同时比较大小和斜率，默认等权
options_extra.weight_list = [1,1];

[FITobj, SubsIndexL] = fitprepare(hr_nsoc, H_soc_sym);
Loss_func_TB = @(para) TBkit.loss_func(para, ...
    'FITobj','FITobj',...
    'DFTBAND','EIG_DFT_soc',...
    'extra','options_extra',...
    'algorithm','pure_comparison',...
    'SubsIndexL',SubsIndexL ...
);   

init_guess = [0.1 0.1];

options = optimset('PlotFcns',@optimplotfval,'Display','iter');
res = fminsearch(Loss_func_TB, init_guess, options);
disp(lambda_syms)
disp(res)
%% plot
H_soc_num_udud = subs(H_soc_sym_udud, lambda_syms, res);
hr_soc_fitted = hr_nsoc;
hr_soc_fitted.HnumL(:,:,hr_soc_fitted.Line_000) = hr_soc_fitted.HnumL(:,:,hr_soc_fitted.Line_000) + H_soc_num_udud;

EIG_hr_fitted = hr_soc_fitted.EIGENCAR_gen();

bandplot({EIG_DFT_soc, EIG_hr_fitted}, [-3,3], 'Color', [1 0 0; 0 0 1]);
legend(["DFT", "wannier"])
title("after fitting")
%% udud to uudd and write to file
if 0 == 1
    WAN_NUM = hr_nsoc.WAN_NUM;
    udud2uudd = [1:2:(WAN_NUM-1),2:2:(WAN_NUM)];
    H_soc_num_uudd = H_soc_num_udud(udud2uudd,udud2uudd);
    
    fid = fopen('soc_term_on_site.dat', 'w');
    for k = 1:numel(H_soc_num_uudd)
        fprintf(fid, '%.16e %.16e\n', real(H_soc_num_uudd(k)), imag(H_soc_num_uudd(k)));
    end
    fclose(fid);
end
