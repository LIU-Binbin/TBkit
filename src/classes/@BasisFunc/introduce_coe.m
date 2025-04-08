function [BFuncL] = introduce_coe(BFuncL,coeL)
for i =1:numel(BFuncL)
BFuncL(i).coe = BFuncL(i).coe*coeL(i);
end
end
