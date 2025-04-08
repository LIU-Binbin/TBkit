function Am = rotate(A,Rc,Rf,tf,rightorleft,optionsOper,options)
arguments
A BasisFunc;
Rc {Spin.mustBeSize(Rc,[3 3])}= diag([1 1 1]);
Rf {Spin.mustBeSize(Rf,[3 3])}= diag([1 1 1]);
tf {Spin.mustBeSize(tf,[3 1;1 3])}= ([0 0 0]);
rightorleft = 'right';
optionsOper.Oper = [];
options.sym = false;
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
if ~isempty(optionsOper.Oper)
Rc = optionsOper.Oper.R;
Rf = optionsOper.Oper.Rf;
tf = optionsOper.Oper.tf;
options.conjugate = optionsOper.Oper.conjugate;
options.antisymmetry = optionsOper.Oper.antisymmetry;
optionsCell = namedargs2cell(options);
Am = rotate(A,Rc,Rf,tf,rightorleft,optionsCell{:});
return
end
optionsCell = namedargs2cell(options);
Am = rotaterow(A(1,:),Rc,Rf,tf,rightorleft,optionsCell{:});
for i = 2:size(A,1)
Am = [Am;rotaterow(A(i,:),Rc,Rf,tf,rightorleft,optionsCell{:})];
end
end
