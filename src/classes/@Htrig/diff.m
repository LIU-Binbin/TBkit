function H_htrig = diff(H_htrig,dir,options)
arguments
H_htrig Htrig;
dir = 1;
options.Accuracy = 1e-6;
end
switch H_htrig.Type
case 'sincos'
VarUsing = H_htrig.VarsSeqLcart(1:H_htrig.Dim);
K = VarUsing;
HsymL_trig_tmp = diff(H_htrig.HsymL,K(dir));
AdditionCoe = zeros(H_htrig.Kinds,1,'sym');
for i  = 1:H_htrig.Kinds
[AdditionCoeTmp,HsymL_trigTmp] = coeffs(HsymL_trig_tmp(i));
if isempty(AdditionCoeTmp) && isempty(HsymL_trigTmp)
AdditionCoe(i) = 0;
HsymL_trig_tmp(i) = 1;
else
AdditionCoe(i) = fold(@times,AdditionCoeTmp);
HsymL_trig_tmp(i) =  fold(@times,HsymL_trigTmp);
end
end
if H_htrig.num
H_htrig.HsymL_trig = HsymL_trig_tmp;
H_htrig.HnumL = TBkit.matrixtimespage(double(AdditionCoe),H_htrig.HnumL);
else
H_htrig.HsymL_trig = HsymL_trig_tmp;
H_htrig.HcoeL = TBkit.matrixtimespage(AdditionCoe,H_htrig.HcoeL);
end
case 'exp'
syms k_x k_y k_z real;
K = [k_x;k_y;k_z];
AdditionCoe = diff(H_htrig.HsymL,K(dir))./H_htrig.HsymL;
if H_htrig.num
H_htrig.HnumL = TBkit.matrixtimespage(double(AdditionCoe),H_htrig.HnumL);
else
H_htrig.HcoeL = TBkit.matrixtimespage(AdditionCoe,H_htrig.HcoeL);
end
case 'mat'
if H_htrig.num
H_htrig.HnumL = 1i*TBkit.matrixtimespage(H_htrig.HsymL_numL(:,dir),H_htrig.HnumL);
else
H_htrig.HcoeL = 1i*TBkit.matrixtimespage(H_htrig.HsymL_coeL(:,dir),H_htrig.HcoeL);
end
case 'list'
if H_htrig.num
H_htrig.HnumL = 1i*H_htrig.HnumL.*H_htrig.HsymL_numL(:,dir);
else
H_htrig.HcoeL = 1i*H_htrig.HcoeL.*H_htrig.HsymL_coeL(:,dir);
end
end
H_htrig =  simplify(H_htrig,options.Accuracy);
end
