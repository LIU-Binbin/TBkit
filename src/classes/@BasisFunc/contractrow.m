function B = contractrow(A,options)
arguments
A BasisFunc;
options.forgetcoe = false;
options.fast = true;
options.hybird = false;
options.spincoupled = false;
options.orbcoupled = false;
options.raw = true;
options.sym = false;
options.conjugate = false;
options.antisymmetry = false;
end
B = A;
if options.hybird
error('not support yet');
return
end
if iscell(A(1).BFuncL)
error('not support yet');
return;
end
A = cleanrow(A);
if  options.spincoupled
end
if options.orbcoupled
end
if ~options.spincoupled && ~options.orbcoupled
BFuncLtmp = ([A.BFuncL]);
coeLtmp = ([A.coe]);
BFuncLtmp = BasisFunc.introduce_coe(BFuncLtmp,coeLtmp);
BFuncLtmp = contractrow(BFuncLtmp);
[coeLtmp,BFuncLtmp] = BasisFunc.extract_coe(BFuncLtmp);
switch class(BFuncLtmp)
case 'Qnum'
spinL = A(1).spin;
otherwise
spinL = A(1).spin;
end
orbL = A(1).BForb;
B = BasisFunc(BFuncLtmp,spinL,1,coeLtmp,orbL);
end
B = cleanrow(B);
end
