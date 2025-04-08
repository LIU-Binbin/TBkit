function [EIGENCAR,WAVECAR,WEIGHTCAR] = EIGENCAR_gen_wire(H_htrig,options)
arguments
H_htrig Htrig;
options.fermi double = 0;
options.norb double = -1;
options.klist  = H_htrig.klist_cart;
options.para  = [];
options.paraname ;
options.issym  logical =false;
options.spdB double = 0;
end
switch H_htrig.Type
case 'sincos'
fprintf('use normal eigencar gen or ?\n');
return;
otherwise
end
fermi = options.fermi;
norb_enforce  = options.norb;
if isempty(H_htrig.klist_cart) && ~isnumeric(options.klist)
if isstr(options.klist)
H_htrig = H_htrig.kpathgen3D(options.klist);
else
try
H_htrig = H_htrig.kpathgen3D('KPOINTS_wire');
catch
end
end
options.klist = H_htrig.klist_cart;
end
klist_cart_tmp = options.klist;
Hnum_list = H_htrig.HnumL ;
NSLAB = (H_htrig.Nslab ==0) + H_htrig.Nslab;
NS = prod(NSLAB);
NWAVE_origin =  H_htrig.Basis_num* NS;
NWAVE = H_htrig.Basis_num* NS-sum(H_htrig.rm_list);
orb_list  = H_htrig.orbL;
HcoeList = H_htrig.HcoeL ;
VarUsing = H_htrig.VarsSeqLcart(1:H_htrig.Dim);
if isempty(H_htrig.rm_list)
HSVCAR_hinge = TBkit.HSVCAR_gen(orb_list,'hinge',0.05,[0.5,0.5,0.5],3);
else
HSVCAR_hinge = TBkit.HSVCAR_gen(orb_list,'surf',0.05,[0.5,0.5,0.5],3);
end
signlist = sign((orb_list(:,1)-0.5).*(orb_list(:,2)-0.5));
if isempty(options.para)
if H_htrig.Basis_num > 1000
print_mode = 1;
else
print_mode = 0;
end
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
Hmat{numel(H_htrig.HsymL_trig)} = sparse(NWAVE_origin,NWAVE_origin);
for i = 1:numel(H_htrig.HsymL_trig)
try
Hfuntemp = matlabFunction(HcoeList(:,:,i),'Vars',VarUsing);
catch
error('You have not subs the function');
end
Hmat{i} = kron(H_htrig.Hmat_pre{i},Hfuntemp(Input{:}));
end
Hout = fold(@plus,Hmat);
Hout = tril(Hout,-1)+diag(real(diag((Hout))))+tril(Hout,-1)';
if ~isempty(H_htrig.rm_list)
Hout(H_htrig.rm_list,:) = [];
Hout(:,H_htrig.rm_list) = [];
end
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
if options.issym
Hout = Hout+eye()*options.spdB;
end
[A, U]=eigs(Hout,NBANDS,options.spdB+fermi+1e-6,'IsSymmetricDefinite',options.issym);
[A, U]=park.sorteig(U,A);
else
end
EIGENCAR(:,ki) = diag(U);
WAVECAR(:,:,ki) = A;
[~,WEIGHTCAR(:,ki)] = Htrig.COLORCAR_gen(A,HSVCAR_hinge,signlist);
if print_mode ==1
fprintf('%d th kpoints has been caculated in %d kpoints total\n',ki,kn);
end
end
else
Npara = size(options.para ,1);
paraN = size(options.para ,2);
if H_htrig.Basis_num > 500
print_mode = 1;
else
print_mode = 0;
end
[kn,~] = size(klist_cart_tmp);
if kn ==1
tmpmode = 'para';
else
tmpmode = 'wire';
end
if norb_enforce <0
NBANDS=NWAVE;
elseif norb_enforce >0
NBANDS=norb_enforce;
else
end
if strcmp(tmpmode,'wire')
WAVECAR  = [];
EIGENCAR{Npara} = zeros(NBANDS,kn);
WEIGHTCAR{Npara} = zeros(NBANDS,kn);
for j = 1:Npara
fprintf('**************************************************************************************\n');
for i = 1:paraN
fprintf('%s :',mat2str(string(sym(options.paraname(i)))));
fprintf('%f\n',options.para(j,i));
end
% fprintf('**************************************************************************************\n');
EIGENCAR_tmp = zeros(NBANDS,kn);
WEIGHTCAR_tmp = zeros(NBANDS,kn);
if strcmp(H_htrig.Type,'sparse')
for i = 1:H_htrig.Kinds
Hnum_list{i} = subs(H_htrig.HnumL{i},sym(options.paraname),options.para(j,:));
end
else
for i =1:numel( H_htrig.HsymL_trig)
temp_str = ["H_fun_t{i} =@(k_x,k_y,k_z) H_htrig.HcoeL(:,:,i)(k_x,k_y,k_z",string(options.para(j,:))];
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
Htemp = Htemp +Hnum_list{i}*double(H_htrig.HsymL_trig(i));
end
Hout = Htemp;
else
Hmat{numel(H_htrig.HsymL_trig)} = sparse(NWAVE_origin,NWAVE_origin);
for i = 1:numel(H_htrig.HsymL_trig)
try
Hfuntemp = H_fun_t{i};
catch
error('You have not subs the function');
end
Hmat{i} = kron(H_htrig.Hmat_pre{i},Hfuntemp(k_x,k_y,k_z));
end
Hout = Hmat{1};
for i = 2:numel(H_htrig.HsymL_trig)
Hout =Hout +Hmat{i};
end
Hout = tril(Hout,-1)+diag(real(diag((Hout))))+tril(Hout,-1)';
end
if ~isempty(H_htrig.rm_list)
Hout(H_htrig.rm_list,:) = [];
Hout(:,H_htrig.rm_list) = [];
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
if options.issym
Hout = Hout+eye()*options.spdB;
end
[A, U]=eigs(Hout,NBANDS,options.spdB+fermi+1e-6,'IsSymmetricDefinite',options.issym);
[A, U]=park.sorteig(U,A);
catch
[A, U]=eig(Hout,'vector');
U = diag(U(NWAVE/2-NBANDS/2+1:NWAVE/2+NBANDS/2));
end
else
end
EIGENCAR_tmp(:,ki) = diag(U);
[~,WEIGHTCAR_tmp(:,ki)] = Htrig.COLORCAR_gen(A,HSVCAR_hinge,signlist);
if print_mode ==1
fprintf('%d th kpoints has been caculated in %d kpoints total\n',ki,kn);
end
end
EIGENCAR{j} = EIGENCAR_tmp;
WEIGHTCAR{para} = WEIGHTCAR_tmp;
end
else
WAVECAR  = zeros(NWAVE,NBANDS,Npara);
EIGENCAR = zeros(NBANDS,Npara);
WEIGHTCAR = zeros(NBANDS,Npara);
k_x=klist_cart_tmp(1);
k_y=klist_cart_tmp(2);
k_z=klist_cart_tmp(3);
for j = 1:Npara
for i = 1:paraN
fprintf('%s :',mat2str(string(sym(options.paraname(i)))));
fprintf('%f\n',options.para(j,i));
end
TheHcoeL = H_htrig.HcoeL;
TheHcoeL =  subs(TheHcoeL,sym(options.paraname),options.para(j,:));
for i =1:numel( H_htrig.HsymL_trig)
H_fun_t{i} = matlabFunction(TheHcoeL(:,:,i),'Vars',H_htrig.VarsSeqLcart(1:H_htrig.Dim));
end
if strcmp(H_htrig.Type,'sparse')
Htemp=sparse(H_htrig.Basis_num ,H_htrig.Basis_num);
for i=1:H_htrig.Kinds
Htemp = Htemp +Hnum_list{i}*double(H_htrig.HsymL_trig(i));
end
Hout = Htemp;
else
Hmat{numel(H_htrig.HsymL_trig)} = sparse(NWAVE_origin,NWAVE_origin);
for i = 1:numel(H_htrig.HsymL_trig)
try
Hfuntemp = H_fun_t{i};
catch
error('You have not subs the function');
end
Hmat{i} = kron(H_htrig.Hmat_pre{i},Hfuntemp(k_x,k_y,k_z));
end
Hout = Hmat{1};
for i = 2:numel(H_htrig.HsymL_trig)
Hout =Hout +Hmat{i};
end
Hout = tril(Hout,-1)+diag(real(diag((Hout))))+tril(Hout,-1)';
end
if ~isempty(H_htrig.rm_list)
Hout(H_htrig.rm_list,:) = [];
Hout(:,H_htrig.rm_list) = [];
end
if norb_enforce <0
try
[A, U]=eig(full(Hout));
catch
disp([k_x,k_y,k_z] );
disp(Hout);
disp(H_htrig.Hfun);
error('check this k point');
end
elseif norb_enforce >0
if options.issym
Hout = Hout+eye()*options.spdB;
end
[A, U]=eigs(Hout,NBANDS,options.spdB+fermi+1e-6,'IsSymmetricDefinite',options.issym);
[A, U]= park.sorteig(U,A);
else
end
EIGENCAR(:,j) = diag(U);
[~,WEIGHTCAR(:,j)] = Htrig.COLORCAR_gen(A,HSVCAR_hinge,signlist);
WAVECAR(:,:,j) = A;
if print_mode ==1
fprintf('%d th para has been caculated in %d Npara total\n',j,Npara);
end
end
end
end
end
