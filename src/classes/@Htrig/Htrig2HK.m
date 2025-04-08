function H_hk = Htrig2HK(H_htrig,kpoints_f,options)
arguments
H_htrig Htrig;
kpoints_f = [0,0,0];
options.sym =false;
options.Order = 1;
end
syms k_x k_y k_z real;
H_hk = HK(H_htrig.Basis_num,options.Order);
H_hk = H_hk.TBkitCopy(H_htrig);
if  options.sym
H_hk.Rm = sym(H_hk.Rm);
kpoints_r = kpoints_f * H_hk.Gk;
fprintf('please check whether the kpoint(cartesian) is properly set.');
disp(kpoints_r);
else
kpoints_r = kpoints_f * H_hk.Gk;
end
pb = TBkit_tool_outer.CmdLineProgressBar('Transforming ');
for i = 1:H_htrig.Kinds
pb.print(i,H_htrig.Kinds,' th HsymL_trig into HK obj ...');
tmp_HsymL_trig = subs(H_htrig.HsymL_trig(i),...
[k_x k_y k_z],...
[k_x-kpoints_r(1),k_y-kpoints_r(2),k_z-kpoints_r(3)]);
symbolic_polynomial = taylor(tmp_HsymL_trig,[k_x k_y k_z],'Order',options.Order+1);
H_hk = H_hk.setup_rough(symbolic_polynomial,H_htrig.HcoeL(:,:,i),true);
end
pb.delete();
end
