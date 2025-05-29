%% initial check

EIG_DFT_soc = EIGENVAL_read();
Ef_DFT_soc = 9.9497;
EIG_DFT_soc = EIG_DFT_soc - Ef_DFT_soc;

hr_nsoc = HR.from_wannier90();
EIG_hr = hr_nsoc.EIGENCAR_gen();
Ef_DFT_nsoc = 9.9930;
EIG_hr = EIG_hr(1:40,:) - Ef_DFT_nsoc;


bandplot( {EIG_DFT_soc,EIG_hr}, 'Color', [1 0 0; 0 0 1]);
legend(["DFT soc", "hr nsoc"])
%% wannier info

element_names = ["Mn", "Pt"];
element_atom_nums = [2 2];
element_projs = {2, [0,1,2]};

[H_soc_full, lambda_syms] = soc_term_udud_add( ...
    'element_names', element_names, ...
    'element_atom_nums',element_atom_nums, ...
    'element_projs',element_projs);
[orbL, quantumL, elementL] = wout_read
% from wannier90.wout & Rm;

% elementL = [ones(2*2*5,1)*25;ones(2*2*(1+3+5),1)*78];
% qnumL(:,1) = [ones(2*2*5,1)*4;ones(2*2*(1+3+5),1)*6];
% qnumL(:,4) = kron(ones(1*(2*5 + 2*(1+3+5)),1),[1/2;-1/2]);
% qnumL(:,2) = [
%     ones(1*(2*2*5),1)*2;
%     ones(1*(2*1),1)*0;
%     ones(1*(2*3),1)*1;
%     ones(1*(2*5),1)*2;
%     ones(1*(2*1),1)*0;
%     ones(1*(2*3),1)*1;
%     ones(1*(2*5),1)*2;
%     ];
% qnumL(:,3) = [
%     0;0;1;1;-1;-1;2;2;-2;-2;
%     0;0;1;1;-1;-1;2;2;-2;-2;
%     0;0;
%     0;0;1;1;-1;-1;
%     0;0;1;1;-1;-1;2;2;-2;-2;
%     0;0;
%     0;0;1;1;-1;-1;
%     0;0;1;1;-1;-1;2;2;-2;-2;
%     ];
[H_soc_full2, lambda_syms2] = soc_term_udud_add(elementL, quantumL, 'mode','direct');
%%
H_soc_full2
hr_soc = hr_nsoc;

% hr_soc.HcoeL = sym(hr_soc.HnumL);
tic
%H_soc_full2 = subs(H_soc_full2,lambda_syms2(3),0);
%初值给定
toc

%
tic
[FITobj,SubsIndexL] = fitprepare(hr_nsoc,H_soc_full2);

toc
%% search
%拟合参数设置
%设定拟合范围：
options_extra.NBAND_range_DFT = [13 : (13+40-1)];
options_extra.NBAND_range = [1:40];
options_extra.klist_range = ':';
% options_extra.E_range = [-2,1.5] - Ef_DFT_nsoc;
options_extra.E_range = [];
%默认大小和斜率等权
options_extra.weight_list = [1,1];
%键值对：
%  ‘EIGENCAR_DFT’ 为你的DFT 能带变量名
% FITobj 为你的拟合对象变量名
%  ‘extra’ 为你的设定的能量范围
%  'algorithm' 为你选择的比较方法 默认为同时比较大小和斜率
 

Loss_func_TB = @(para) TBkit.loss_func(para, ...
    'FITobj','FITobj',...
    'DFTBAND','EIG_DFT_soc',...
    'extra','options_extra',...
    'algorithm','pure_comparison',...
    'SubsIndexL',SubsIndexL ...
)          

FITobj.symvar_list
x0 = [0.1 0.4,0.01];
%进行拟合
options = optimset('PlotFcns',@optimplotfval,'Display','iter');
x = fminsearch(Loss_func_TB,x0,options);
% x = fmincon(Loss_func_TB,x0,...
%     [0 0 1],0.1,[],[],[],[],[],...
%    options);
[FITobj_n,EQ] = FITobj.subs(x,'SubsIndex',SubsIndexL)
% FITobj_n.HnumL(SubsIndexL) = double(subs(FITobj.HcoeL,FITobj.symvar_list,x));
% EQL =(FITobj.symvar_list == x)

EIG_hr = FITobj_n.EIGENCAR_gen()-Ef_DFT_nsoc;
% 画出能带
bandplot( {EIG_DFT_soc,EIG_hr},[-1,1], 'Color', [1 0 0; 0 0 1]);
legend(["DFT soc", "hr nsoc"])
return;
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
