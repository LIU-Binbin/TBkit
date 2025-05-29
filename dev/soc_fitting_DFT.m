%% initial check
EIG_DFT_soc = EIGENVAL_read();
Ef_DFT_soc = 9.9497;
EIG_DFT_soc = EIG_DFT_soc(13:12+40,:) - Ef_DFT_soc;

hr_nsoc = HR.from_wannier90();
EIG_hr = hr_nsoc.EIGENCAR_gen();
Ef_DFT_nsoc = 9.9930;
EIG_hr = EIG_hr(1:40,:) - Ef_DFT_nsoc;
%% 高能带拟合不太好，截取下面的能带
bandplot( {EIG_DFT_soc,EIG_hr}, 'Color', [1 0 0; 0 0 1]);
legend(["DFT soc", "hr nsoc"])
%% wannier info
element_names = ["Mn", "Pt"];
element_atom_nums = [2 2];
element_projs = {2, [0,1,2]};

[H_soc_full, lambda_syms] = soc_term_udud_add(element_names, element_atom_nums, element_projs);
%% search
soc_fitting_handle = @(lambda_nums) soc_fitting(lambda_nums, lambda_syms, EIG_DFT_soc, hr_nsoc, H_soc_full);
lambda_guess = [0.2 0.2];
[lambda_out, loss_out]= fminsearch(soc_fitting_handle, lambda_guess);
%% check
for i = 1:length(lambda_out)
    disp(string(lambda_syms(i))+" = "+lambda_out(i)+" eV")
end

H_soc_num = subs(H_soc_full, lambda_syms, lambda_out);
hr_nsoc.HnumL(:,:,hr_nsoc.Line_000) = hr_nsoc.HnumL(:,:,hr_nsoc.Line_000) + H_soc_num;
EIG_hr = hr_nsoc.EIGENCAR_gen();
Ef_DFT_nsoc = 9.9930;
EIG_hr = EIG_hr(1:40,:) - Ef_DFT_nsoc;
bandplot( {EIG_DFT_soc,EIG_hr}, 'Color', [1 0 0; 0 0 1]);
%% udud to uudd and write to file
WAN_NUM = hr_nsoc.WAN_NUM;
udud2uudd = [1:2:(WAN_NUM-1),2:2:(WAN_NUM)];
H_soc_num_uudd = H_soc_num(udud2uudd,udud2uudd);

fid = fopen('soc_term_on_site.dat', 'w');
for k = 1:numel(H_soc_num_uudd)
    fprintf(fid, '%.16e %.16e\n', real(H_soc_num_uudd(k)), imag(H_soc_num_uudd(k)));
end
fclose(fid);
%%
function loss = soc_fitting(params, params_sym, EIG_DFT_soc, hr_nsoc, H_soc_sym)
H_soc_num = subs(H_soc_sym, params_sym, params);

hr_nsoc.HnumL(:,:,hr_nsoc.Line_000) = hr_nsoc.HnumL(:,:,hr_nsoc.Line_000) + H_soc_num;

EIG_hr = hr_nsoc.EIGENCAR_gen();
Ef_DFT_nsoc = 9.9930;
EIG_hr = EIG_hr(1:40,:) - Ef_DFT_nsoc;

loss = sum( abs(EIG_DFT_soc - EIG_hr), "all");
end
