function [EIGENCAR,WAVECAR,WEIGHTCAR] = EIGENCAR_gen_slab(H_htrig,options)
arguments
H_htrig Htrig;
options.fermi double = 0;
options.norb double = -1;
options.klist  = H_htrig.klist_cart;
options.para  = [];
options.paraname ;
end
switch H_htrig.Type
case 'sincos'
fprintf('use normal eigencar gen or ?\n');
return;
otherwise
end
fermi = options.fermi;
norb_enforce  = options.norb;
if isempty(H_htrig.klist_cart)
if isstr(options.klist)
H_htrig = H_htrig.kpathgen3D(options.klist);
else
H_htrig = H_htrig.kpathgen3D('KPOINTS_slab');
end
options.klist = H_htrig.klist_cart;
end
klist_cart_tmp = options.klist;
HcoeList = H_htrig.HcoeL ;
VarUsing = H_htrig.VarsSeqLcart(1:H_htrig.Dim);
NSLAB = (H_htrig.Nslab ==0) + H_htrig.Nslab;
NS = NSLAB(1)* NSLAB(2)* NSLAB(3);
NWAVE = H_htrig.Basis_num* NS;
orb_list = H_htrig.orbL;
if H_htrig.Nslab(1) > 1
signlist = sign((orb_list(:,1)-0.5));
dir = 1;
elseif H_htrig.Nslab(2) > 1
dir = 2;
signlist = sign((orb_list(:,2)-0.5));
elseif H_htrig.Nslab(3) > 1
dir = 3;
signlist = sign((orb_list(:,3)-0.5));
elseif H_htrig.Nslab(4) > 1
dir = 4;
signlist = sign((orb_list(:,4)-0.5));
end
HSVCAR_slab = TBkit.HSVCAR_gen(orb_list,'slab',0.05,[0.5,0.5,0.5],dir);
signlist(HSVCAR_slab(:,1) == 0) = 0;
if size(H_htrig.orbL,1) > 1000
print_mode = true;
else
print_mode = false;
end
if isempty(options.para)
[kn,~] = size(klist_cart_tmp);
if norb_enforce <0
NBANDS=NWAVE;
elseif norb_enforce >0
NBANDS=norb_enforce;
else
end
WAVECAR  = zeros(NWAVE,NBANDS,kn);
EIGENCAR = zeros(NBANDS,kn);
WEIGHTCAR = zeros(NBANDS,kn);
for ki =1:kn
k_d=klist_cart_tmp(ki,:);
Input = num2cell(k_d);
Hmat = zeros(NWAVE,NWAVE,numel(H_htrig.HsymL_trig));
for i = 1:numel(H_htrig.HsymL_trig)
try
Hfuntemp = matlabFunction(HcoeList(:,:,i),'Vars',VarUsing);
catch
error('You have not subs the function');
end
Hmat(:,:,i) = kron(H_htrig.Hmat_pre{i},Hfuntemp(Input{:}));
end
Hout = sum(Hmat,3);
Hout = tril(Hout,-1)+diag(real(diag((Hout))))+tril(Hout,-1)';
if norb_enforce <0
try
[A, U]=eig(full(Hout));
catch
disp([ki,Input] );
disp(Hout);
disp(H_htrig.Hfun);
error('check this k point');
end
elseif norb_enforce >0
[A, U]=eigs(Hout,NBANDS,fermi);
[A, U]=park.sorteig(U,A);
else
end
EIGENCAR(:,ki) = diag(U);
WAVECAR(:,:,ki) = A;
[~,WEIGHTCAR(:,ki)] = Htrig.COLORCAR_gen(A,HSVCAR_slab,signlist);
if print_mode ==1
fprintf('%d th kpoints has been caculated in %d kpoints total\n',ki,kn);
end
end
else
Npara = size(options.para ,1);
paraN = size(options.para ,2);
[kn,~] = size(klist_cart_tmp);
if norb_enforce <0
NBANDS=NWAVE;
elseif norb_enforce >0
NBANDS=norb_enforce;
else
end
WEIGHTCAR = [];
WAVECAR  = [];
EIGENCAR{Npara} = zeros(NBANDS,kn);
for j = 1:Npara
fprintf('**************************************************************************************\n');
for i = 1:paraN
fprintf('%s :',mat2str(string(sym(options.paraname(i)))));
fprintf('%f\n',options.para(j,i));
end
% fprintf('**************************************************************************************\n');
EIGENCAR_tmp = zeros(NBANDS,kn);
if strcmp(H_htrig.Type,'sparse')
for i = 1:H_htrig.Kinds
HcoeList{i} = subs(H_htrig.HnumL{i},sym(options.paraname),options.para(j,:));
end
else
for i =1:numel( H_htrig.HsymL_trig)
temp_str = ["H_fun_t{i} =@(k_x,k_y,k_z) H_htrig.Hfun{i}(k_x,k_y,k_z",string(options.para(j,:))];
temp_str = strjoin(temp_str,',');
temp_str = temp_str+");";
eval(temp_str);
end
end
for ki =1:kn
k_x=klist_cart_tmp(ki,1);
k_y=klist_cart_tmp(ki,2);
k_z=klist_cart_tmp(ki,3);
if strcmp(H_htrig.Type,'sparse')
Htemp=sparse(H_htrig.Basis_num ,H_htrig.Basis_num);
for i=1:H_htrig.Kinds
Htemp = Htemp +HcoeList{i}*double(H_htrig.HsymL_trig(i));
end
Hout = Htemp;
else
Hmat = zeros(NWAVE,NWAVE,numel(H_htrig.HsymL_trig));
for i = 1:numel(H_htrig.HsymL_trig)
try
Hfuntemp = H_fun_t{i};
catch
error('You have not subs the function');
end
Hmat(:,:,i) = kron(H_htrig.Hmat_pre{i},Hfuntemp(k_x,k_y,k_z));
end
Hout = sum(Hmat,3);
Hout = tril(Hout,-1)+diag(real(diag((Hout))))+tril(Hout,-1)';
end
if norb_enforce <0
try
[A, U]=eig(full(Hout));
catch
disp([ki,k_x,k_y,k_z] );
disp(Hout);
disp(H_htrig.Hfun);
error('check this k point');
end
elseif norb_enforce >0
try
[A, U]=eigs(Hout,NBANDS,fermi);
[A, U]=park.sorteig(U,A);
catch
[A, U]=eig(Hout,'vector');
U = diag(U(NWAVE/2-NBANDS/2+1:NWAVE/2+NBANDS/2));
end
else
end
EIGENCAR_tmp(:,ki) = diag(U);
if print_mode ==1
fprintf('%d th kpoints has been caculated in %d kpoints total\n',ki,kn);
end
end
EIGENCAR{j} = EIGENCAR_tmp;
end
end
end
