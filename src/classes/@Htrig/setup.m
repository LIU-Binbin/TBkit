function H_htrig = setup(H_htrig,Var_cell,k_cell,mat_cell,silence)
if nargin < 5
silence = false;
end
BASIS_NUM = H_htrig.Basis_num;
if length(Var_cell)~=length(k_cell) && length(k_cell)~=length(mat_cell)
error('error!');
end
if ~silence
pb = TBkit_tool_outer.CmdLineProgressBar('Setting ');
end
nVar_cell = length(Var_cell);
for i =1:nVar_cell
Var = Var_cell{i};
k_symbol = (k_cell{i});
matcell = mat_cell{i};
Kind = H_htrig.k_symbol2Kind(k_symbol);
if isempty(Kind)
Kind = H_htrig.Kinds+1;
H_htrig.HsymL_trig(Kind) = k_symbol;
H_htrig.HcoeL(:,:,Kind) = sym(zeros(BASIS_NUM,BASIS_NUM,1));
H_htrig.HnumL(:,:,Kind)  = (zeros(BASIS_NUM,BASIS_NUM,1));
end
if ~silence
pb.print(i,nVar_cell,'Htrig ...');
end
switch class(Var)
case 'sym'
H_htrig.HcoeL(:,:,Kind) = H_htrig.HcoeL(:,:,Kind)+ matcell*Var;
case 'string'
case 'double'
H_htrig.HnumL(:,:,Kind) = H_htrig.HnumL(:,:,Kind)+ matcell*Var;
end
end
if ~silence
pb.delete();
end
end
