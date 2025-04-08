function [coeL,BFuncL] = extract_coe(BFuncL,options)
arguments
BFuncL
options.sym = false;
options.vpalevel = 6;
end
coeL = ones(size(BFuncL),class(BFuncL(1).coe));
oneterm  = ones(1,1,class(BFuncL(1).coe));
for i =1:numel(BFuncL)
coeL(i) = BFuncL(i).coe;
BFuncL(i).coe =  oneterm;
end
if options.sym
else
coeL = round(coeL,options.vpalevel);
end
end
