function varargout = EIGENCAR_gen(H_hk,options)
arguments
H_hk HK;
options.fermi double = 0;
options.norb double = -1;
options.klist double = H_hk.klist_cart;
options.para  = [];
options.paraname ;
options.show = false;
options.OriginGk = [0 0 0];
options.ax = handle([]);
options.printmode = false;
end
fermi = options.fermi;
norb_enforce  = options.norb;
if isempty(options.klist)
H_hk = H_hk.kpathgen3D('KPOINTS');
klist_cart_tmp = H_hk.klist_cart;
else
klist_cart_tmp = options.klist;
end
klist_cart_tmp = klist_cart_tmp + options.OriginGk;
if options.show
if iisempty(options.ax)
ax = TBkit.BZplot(H_hk.Rm,'color','r');
else
ax = options.ax;
end
end
if H_hk.Basis_num > 500
print_mode = 1;
else
print_mode = 0;
end
[kn,~] = size(klist_cart_tmp);
if norb_enforce <0
NBANDS=H_hk.Basis_num;
elseif norb_enforce >0
NBANDS=norb_enforce;
else
end
WAVECAR  = zeros(H_hk.Basis_num,NBANDS,kn);
EIGENCAR = zeros(NBANDS,kn);
HsymL_fun = matlabFunction( H_hk.HsymL,'Vars',H_hk.VarsSeqLcart(1:H_hk.Dim));
for ki =1:kn
Input = num2cell(klist_cart_tmp(ki,:));
kL = HsymL_fun(Input{:});
Hout = H_hk.HnumL;
for i =1:H_hk.Kinds
Hout(:,:,i) = Hout(:,:,i).*kL(:,i);
end
Hout = sum(Hout,3);
Hout = (Hout+Hout')/2;
if norb_enforce <0
try
[A, U]=eig(full(Hout));
catch
disp([ki,x,y,z] );
disp(Hout);
disp(H_hk.Hfun);
error('check this k point');
end
elseif norb_enforce >0
[A, U]=eigs(Hout,NBANDS,fermi);
[A, U]= park.sorteig(U,A);
else
end
EIGENCAR(:,ki) = diag(U);
WAVECAR(:,:,ki) = A;
if print_mode ==1
fprintf('%d th kpoints has been caculated in %d kpoints total\n',ki,kn);
end
end
varargout{1} = EIGENCAR;
varargout{2} = WAVECAR;
if options.show
[varargout{3}] = TBkit.klist_show(...
'klist',klist_cart_tmp,...
'ax',ax);
end
end
