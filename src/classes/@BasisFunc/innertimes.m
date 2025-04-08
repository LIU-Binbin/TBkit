function C = innertimes(A,B,options)
arguments
A BasisFunc;
B BasisFunc;
options.forgetcoe = false;
options.fast = true;
options.hybird = false;
options.spincoupled = false;
options.orbcoupled = false;
options.raw = true
end
optionsCell = namedargs2cell(options);
if ~options.hybird && ~options.spincoupled && ~options.orbcoupled
coeLtmp0 = A.coe*B.coe;
BFuncLtmp = (A.BFuncL .* B.BFuncL);
[coeLtmp,BFuncLtmp] = BasisFunc.extract_coe(BFuncLtmp);
spinL = A(1).spin;
orbL = A(1).BForb;
C = BasisFunc(BFuncLtmp,spinL,1,coeLtmp0*coeLtmp,'orbL',orbL);
end
end
