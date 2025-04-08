function A_Lj = rotaterow(A,Rc,Rf,tf,rightorleft,options)
arguments
A BasisFunc;
Rc {Spin.mustBeSize(Rc,[3 3])}= diag([1 1 1]);
Rf {Spin.mustBeSize(Rf,[3 3])}= diag([1 1 1]);
tf {Spin.mustBeSize(tf,[3 1;1 3])}= ([0 0 0]);
rightorleft = 'right';
options.sym = true;
options.conjugate = false;
options.antisymmetry = false;
options.forgetcoe = false;
options.fast = true;
options.hybird = false;
options.spincoupled = false;
options.orbcoupled = false;
options.raw = true;
options.vpalevel = 6;
options.center = [0,0,0];
end
optionsCell = namedargs2cell(options);
if options.hybird
error('not support yet');
return
else
BFuncLtmp = ([A.BFuncL]);
coeLtmp = ([A.coe]);
BFuncLtmp = BasisFunc.introduce_coe(BFuncLtmp,coeLtmp);
end
if iscell(A(1).BFuncL)
error('not support yet');
return;
end
if ~options.orbcoupled
orbL = BasisFunc.rotation_orb(A(1).BForb,Rf.',tf,optionsCell{:});
if ~options.spincoupled
if isa(A(1).BFuncL,'Qnum')
BFuncLtmp = rotaterow(BFuncLtmp,Rc,rightorleft,...
'sym',options.sym,'antisymmetry',options.antisymmetry,'conjugate',options.conjugate);
BFuncLtmp = contractrow(BFuncLtmp);
spinL = Spin([BFuncLtmp.s],[BFuncLtmp.sz]);
[coeLtmp,BFuncLtmp] = BasisFunc.extract_coe(BFuncLtmp,'sym',options.sym,'vpalevel',options.vpalevel);
A_Lj = BasisFunc(BFuncLtmp,spinL,1,coeLtmp,orbL);
else
BFuncLtmp = rotaterow(BFuncLtmp,Rc,rightorleft,...
'sym',options.sym,'antisymmetry',options.antisymmetry,'conjugate',options.conjugate);
spinL = rotaterow([A.spin],Rc,rightorleft,...
'sym',options.sym,'antisymmetry',options.antisymmetry,'conjugate',options.conjugate);
BFuncLtmp = contractrow(BFuncLtmp);
spinL = contractrow(spinL);
end
else
end
else
error('not imply yet');
end
end
