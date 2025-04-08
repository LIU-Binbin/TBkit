function H_htrig = add_empty_one(H_htrig,vector)
if isempty(H_htrig.coe)||isempty(H_htrig.num)
[~,~,H_htrig] = H_htrig.NumOrCoe();
end
switch H_htrig.Type
case {'mat','list'}
if H_htrig.num
vector = double(vector);
for i = 1:size(vector,1)
vector_single = (vector(i,:));
try
if (ismember(vector_single,H_htrig.HsymL_numL,'rows'))
continue;
end
catch
end
Kind = H_htrig.Kinds +1;
H_htrig.HsymL_numL(Kind,:) = (vector_single);
H_htrig.HnumL(:,:,Kind) = zeros(H_htrig.Basis_num);
end
end
if H_htrig.coe
vector = sym(vector);
for i = 1:size(vector,1)
vector_single = (vector(i,:));
try
if (ismember(vector_single,H_htrig.HsymL_coeL,'rows'))
continue;
end
catch
end
Kind = H_htrig.Kinds +1;
H_htrig.HsymL_coeL(Kind,:) = (vector_single);
H_htrig.HcoeL(:,:,Kind) = zeros(H_htrig.Basis_num,'sym');
end
end
otherwise
Kind = H_htrig.Kinds+1;
BASIS_NUM = H_htrig.Basis_num;
H_htrig.HsymL_trig(Kind) = vector;
H_htrig.HcoeL(:,:,Kind) = zeros(BASIS_NUM,BASIS_NUM,'sym');
H_htrig.HnumL(:,:,Kind)  = (zeros(BASIS_NUM,BASIS_NUM));
end
end
