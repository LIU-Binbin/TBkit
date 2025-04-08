function H_hk = HR2HK(H_hr,kpoints_frac,options)
arguments
H_hr HR;
kpoints_frac = [0,0,0];
options.sym =false;
options.Order = 2;
options.convention = 'I';
end
if strcmp(H_hr.Type,'list')
H_hr = H_hr.rewrite('rewind',true);
H_hk = HR2HK(H_hr,kpoints_frac,...
'sym',options.sym,...
'Order',options.Order,...
'convention',options.convention);
else
VarsUsing = H_hr.VarsSeqLcart(1:H_hr.Dim);
H_hk = HK(H_hr.WAN_NUM,options.Order,'Dim',H_hr.Dim);
H_hk = H_hk.TBkitCopy(H_hr);
if  options.sym
H_hr.Rm = sym(H_hr.Rm);
H_hr.orbL = sym(H_hr.orbL);
fprintf('please check whether the Rm and orbL are properly set.');
disp(H_hr.Rm);
disp(H_hr.orbL);
end
kpoints_r = kpoints_frac *H_hr.Gk;
vectorList = double(H_hr.vectorL);
if strcmp(options.convention,'II')
for n = 1:H_hr.NRPTS
kpoints_phase = exp(1i*2*pi*(vectorList(n,:))*kpoints_frac.');
tmp_HsymL_trig = exp(1i*(vectorList(n,:))*H_hk.Rm*VarsUsing.');
symbolic_polynomial = taylor(tmp_HsymL_trig,VarsUsing,'Order',options.Order);
H_hk = H_hk.setup_rough(symbolic_polynomial,H_hr.HcoeL(:,:,n)*kpoints_phase);
end
elseif strcmp(options.convention,'I')
H_hr = H_hr.tjmti_gen('sym');
if options.sym
tji_mat =  H_hr.tjmti{1};
else
tji_mat = double(H_hr.tjmti{1});
end
for n = 1:H_hr.NRPTS
for i = 1:H_hr.WAN_NUM
for j = 1:H_hr.WAN_NUM
tjitmp = reshape(tji_mat(i,j,:),[1 H_hr.Dim]);
R_add_t_vector = (vectorList(n,1:H_hr.Dim)*H_hr.Rm+tjitmp).';
kpoints_phase = exp(1i*kpoints_r*(R_add_t_vector));
tmp_HsymL_trig = exp(1i*VarsUsing*(R_add_t_vector));
symbolic_polynomial = taylor(tmp_HsymL_trig,VarsUsing,'Order',options.Order);
H_hk = H_hk.setup_single(simplify(symbolic_polynomial*H_hr.HcoeL(i,j,n)*kpoints_phase),i,j);
end
end
end
else
H_htrig = HR2Htrig(H_hr,'sym',options.sym);
H_hk = H_htrig.Htrig2HK(kpoints_frac,'sym',options.sym,'Order',options.Order);
end
end
H_hk.HcoeL = simplify(H_hk.HcoeL);
end
