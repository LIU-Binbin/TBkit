function H_hr = Htrig2HR(H_htrig,options)
arguments
H_htrig Htrig;
options.POSCAR = 'POSCAR';
end
switch H_htrig.Type
case 'sincos'
Htrig_exp = H_htrig.rewrite();
case 'exp'
Htrig_exp = H_htrig;
otherwise
end
Hexp = Htrig_exp.HcoeL;
WAN_NUM = Htrig_exp.Basis_num;
H_hr = HR(WAN_NUM,'Dim',Htrig_exp.Dim);
H_hr = H_hr.TBkitCopy(Htrig_exp);
if isempty(H_hr.orbL)
warning(['use ',options.POSCAR,'enforcely, please check it']);
H_hr = H_hr < options.POSCAR;
end
H_hr = H_hr.tjmti_gen();
tji_mat_r = H_hr.tjmti{1};
VarsUsing = H_htrig.VarsSeqLcart(1:H_htrig.Dim);
hsym = Htrig_exp.HsymL;
RM = H_hr.Rm;
DIM = H_hr.Dim;
for n = 1:length(hsym)
if isequal(hsym(n), sym(1))
kd_num = zeros(1,DIM);
else
[ChirdrenCell,Type] = TBkit.fixedchildren(hsym(n),'exp_inner');
if strcmp(Type,'sum') || strcmp(Type,'prod')
error('? something wrong? You should Debug here!');
elseif strcmp(Type,'inner')
kd = ChirdrenCell{1}/1i;
end
for d = 1:DIM
cDim{d} = subs(kd,VarsUsing(d),1) - subs(kd,VarsUsing(d),0);
end
kd_num = double(fold(@horzcat,cDim));
end
for i = 1:WAN_NUM
for j = 1:WAN_NUM
if isequal(Hexp(i,j,n), sym(0))
continue
end
kji_num = reshape(tji_mat_r(i,j,:),[1,DIM]);
kr =   kd_num-kji_num;
vector = round(kr/RM);
H_hr = H_hr.set_hop((Hexp(i,j,n)),i,j,vector,'symadd');
end
end
end
end
