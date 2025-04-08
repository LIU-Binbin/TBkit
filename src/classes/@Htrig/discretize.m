function H_htrig2 = discretize(H_htrig,Nslab,options)
arguments
H_htrig Htrig;
Nslab double = [0,10,0];
options.rmfunc function_handle=@()(1);
options.Rotation  = sym(eye(3));
end
if strcmp(functions(options.rmfunc).function , '@()(1)')
rm_mode =false;
else
rm_mode =true;
end
if isequal(options.Rotation ,sym(eye(3)))
rotate_mode = false;
else
rotate_mode = true;
end
if rotate_mode
H_htrig = H_htrig.rotation(options.Rotation);
end
H_htrig2 = Htrig(H_htrig.Basis_num,'Dim',H_htrig.Dim);
H_htrig2.seeds = ["Sigma_x","Sigma_y","Sigma_z","Sigma_w"];
H_htrig2.Type = 'slab';
H_htrig2.Nslab = Nslab;
NSLAB = (H_htrig2.Nslab ==0) + H_htrig2.Nslab;
NS = fold(@times,NSLAB);
case_d = Nslab>1;
seeds = ["x","y","z","w"];
StrSinUsing = "sin(k_" + seeds + ")";
StrCosUsing = "cos(k_" + seeds + ")";
StrSigma_1__N  = "Sigma_" + seeds + "_1__N";
StrSigma_N__N  = "Sigma_" + seeds + "_N__N";
pat_d_pre   = "Sigma_" + seeds + "_";
SymSinUsing = str2sym(StrSinUsing);
SymCosUsing = str2sym(StrCosUsing);
SymSigma_1__N = str2sym(StrSigma_1__N);
SymSigma_N__N = str2sym(StrSigma_N__N);
for d = 1:numel(case_d)
if case_d(d)
H_htrig.HsymL_trig = subs(H_htrig.HsymL_trig,SymSinUsing(d),1i/2 * SymSigma_1__N(d));
H_htrig.HsymL_trig = subs(H_htrig.HsymL_trig,SymCosUsing(d),1 /2 * SymSigma_1__N(d));
pat_d = pat_d_pre(d)+digitsPattern(1)+"__N";
for i =1:length(H_htrig.HsymL_trig)
if ~contains(string(H_htrig.HsymL_trig(i)),pat_d)
H_htrig.HsymL_trig(i) = H_htrig.HsymL_trig(i)*SymSigma_N__N(d);
end
end
H_htrig.Sigmas =[H_htrig.Sigmas;[SymSigma_N__N(d),SymSigma_1__N(d)]];
end
end
H_htrig2.HsymL_trig_bk = [SymSigma_1__N SymSigma_N__N];
H_htrig2.HsymL_trig = sym([]);
count = 0;
for i = 1:H_htrig.Kinds
[coeff_trig,symvar_list_trig,H_htrig2] = split_sym_eq(H_htrig2,H_htrig.HsymL_trig(i));
for j =1:numel(coeff_trig)
count = count+1;
k_cell{count} = symvar_list_trig(j);
mat_cell{count} = H_htrig.HcoeL(:,:,i);
Var_cell{count} = coeff_trig(j);
end
end
H_htrig2 = H_htrig2.setup(Var_cell,k_cell,mat_cell);
orb_tmp = zeros(NS*H_htrig.Basis_num,H_htrig.Dim);
NWAVE  = NS*H_htrig.Basis_num;
if isempty(H_htrig.orbL)
H_htrig.orbL = zeros(H_htrig.Basis_num,H_htrig.Dim);
end
switch H_htrig.Dim
case 1
[i1L            ] = ind2sub(NSLAB,1:NS);
OrbAddL = i1L.' -1;
case 2
[i1L,i2L        ] = ind2sub(NSLAB,1:NS);
OrbAddL = [i1L.' i2L.'] -1;
case 3
[i1L,i2L,i3L    ] = ind2sub(NSLAB,1:NS);
OrbAddL = [i1L.' i2L.' i3L.'] -1;
case 4
[i1L,i2L,i3L,i4L] = ind2sub(NSLAB,1:NS);
OrbAddL = [i1L.' i2L.' i3L.' i4L.'] -1;
end
for i=1:NS
orb_tmp((i-1)*H_htrig.Basis_num+1:i*H_htrig.Basis_num,:) = H_htrig.orbL+OrbAddL(i,:);
end
orb_tmp = orb_tmp./NSLAB;
SizeOrb_tmp = size(orb_tmp);
if rm_mode
try
for d = 1:SizeOrb_tmp(2)
Input{d} = orb_tmp(:,d);
end
H_htrig2.rm_list = options.rmfunc(Input{:});
catch
H_htrig2.rm_list = false(1,NWAVE);
for i =1:NWAVE
for d = 1:SizeOrb_tmp(2)
Input{d} = orb_tmp(i,d);
end
H_htrig2.rm_list(i) = options.rmfunc(Input{:});
end
end
else
end
orb_tmp( H_htrig2.rm_list,:) = [];
H_htrig2.orbL = orb_tmp;
H_htrig2.Hmat_pre{numel(H_htrig2.HsymL_trig)} = sparse(NS,NS);
for i = 1:numel(H_htrig2.HsymL_trig)
H_htrig2.Hmat_pre{i} = H_htrig2.HsymL_trig2mat(H_htrig2.HsymL_trig(i));
end
end
