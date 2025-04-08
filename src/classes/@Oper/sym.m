function [SymMat,SymR] = sym(SymOper)
if numel(SymOper)>1
for i = 1:numel(SymOper)
[SymOper(i).U,SymOper(i).R] = sym(SymOper(i));
end
SymMat = SymOper;SymR = [];
else
SymMat = sym(SymOper.U);
SymR = sym(SymOper.R);
end
end
