function [HoutL] = HCAR_gen(H_htrig,klist,options)
arguments
H_htrig;
klist;
options.Hermi = true;
end
kn = size(klist,1);
switch H_htrig.Type
case 'list'
if isempty(H_htrig.Sparse_vector) && isempty(H_htrig.CutList)
H_htrig = SliceGen(H_htrig);
end
HoutL = zeros(H_htrig.Basis_num,H_htrig.Basis_num,kn);
Hnum_list_k = H_htrig.HnumL.*exp(1i*H_htrig.HsymL_numL(:,1:3)*klist.');
for ki = 1:kn
Hout = zeros(H_htrig.Basis_num,H_htrig.Basis_num);
Hnum_list_ktmp = Hnum_list_k(:,ki);
for i=1:H_htrig.N_Sparse_vector
Hout(H_htrig.Sparse_vector(i,1),H_htrig.Sparse_vector(i,2)) = sum(Hnum_list_ktmp(H_htrig.CutList(i,1):H_htrig.CutList(i,2)));
end
if options.Hermi
HoutL(:,:,ki) = (Hout+Hout')/2;
else
HoutL(:,:,ki) = HoutL;
end
end
case 'sparse'
for ki =1:kn
Htemp=sparse(H_htrig.Basis_num ,H_htrig.Basis_num);
for ki=1:H_htrig.Kinds
Htemp = Htemp +Hnum_list{ki}*double(H_htrig.HsymL_trig(ki));
end
HoutL{ki} = Htemp;
end
case 'mat'
HoutL = zeros(H_htrig.Basis_num,H_htrig.Basis_num,kn);
for ki =1:kn
Factorlist = exp(1i*H_htrig.HsymL_numL*klist(ki,:).');
Hout = sum(TBkit.matrixtimespage(Factorlist,H_htrig.HnumL),3);
if options.Hermi
HoutL(:,:,ki) = (Hout+Hout')/2;
else
HoutL(:,:,ki) = HoutL;
end
end
otherwise
H_htrig.Htrig_num = subs(H_htrig.Htrig_sym);
H_htrig.Hfun = matlabFunction(H_htrig.Htrig_num,'Vars',[sym('k_x'),sym('k_y'),sym('k_z')]);
HoutL = zeros(H_htrig.Basis_num,H_htrig.Basis_num,kn);
for ki =1:kn
Hout = H_htrig.Hfun(klist(ki,1),klist(ki,2),klist(ki,3));
if options.Hermi
HoutL(:,:,ki) = (Hout+Hout')/2;
else
HoutL(:,:,ki) = HoutL;
end
end
end
end
