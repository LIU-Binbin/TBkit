function H_htrig = dualize(H_htrig)
if strcmp(H_htrig.Type,'exp')
BASIS_NUM = H_htrig.Basis_num;
for i=1:length(H_htrig.HsymL_trig)
k_symbol = conj(H_htrig.HsymL_trig(i));
Kind = H_htrig.k_symbol2Kind(k_symbol);
if isempty(Kind)
Kind = H_htrig.Kinds+1;
H_htrig.HsymL_trig(Kind) = k_symbol;
H_htrig.HcoeL(:,:,Kind) = sym(zeros(BASIS_NUM,BASIS_NUM,1));
H_htrig.HnumL(:,:,Kind) = (zeros(BASIS_NUM,BASIS_NUM,1));
end
end
[~,H_htrig.Duality_vector_dist] = ismember(conj(H_htrig.HsymL_trig),H_htrig.HsymL_trig);
end
end
