function H_htrig = simplify(H_htrig,Accuracy,options)
arguments
H_htrig Htrig;
Accuracy = 1e-6;
options.reduce = false;
end
if options.reduce
for i = 1:numel(H_htrig.HcoeL)
if H_htrig.HcoeL(i)~=sym(0)
H_htrig.HcoeL(i) = TBkit.cleanVar(H_htrig.HcoeL,log10(Accuracy));
end
end
end
[~,~,H_htrig] = NumOrCoe(H_htrig);
if H_htrig.coe
H_coeL_tmp = simplify(H_htrig.HcoeL);
H_htrig.HcoeL = H_coeL_tmp;
switch H_htrig.Type
case 'list'
Kinds_list = find(H_coeL_tmp ~=sym(0));
H_htrig = H_htrig.reseq(':',Kinds_list);
case {'mat','exp','sincos'}
zerosMat = zeros(size(H_coeL_tmp(:,:,1)),'sym');
Kinds_list = true(H_htrig.Kinds,1);
for i = 1:H_htrig.Kinds
if isequal(zerosMat,H_coeL_tmp(:,:,i))
Kinds_list(i) = false;
end
end
H_htrig = H_htrig.reseq(':',Kinds_list);
otherwise
end
end
if H_htrig.num
H_numL_tmp = H_htrig.HnumL;
switch H_htrig.Type
case 'list'
Kinds_list = find(abs(H_numL_tmp) > Accuracy);
H_htrig = H_htrig.reseq(':',Kinds_list);
case {'mat','exp','sincos'}
zerosMat = ones(H_htrig.Basis_num)*Accuracy;
Kinds_list = true(H_htrig.Kinds,1);
for i = 1:H_htrig.Kinds
if sum(sum(abs(H_numL_tmp(:,:,i)) > zerosMat))
else
Kinds_list(i) = false;
end
end
H_htrig = H_htrig.reseq(':',Kinds_list);
otherwise
end
end
end
