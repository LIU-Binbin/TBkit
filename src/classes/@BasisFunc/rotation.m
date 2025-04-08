function U = rotation(A,Rc,Rf,tf,optionsConvection,optionsOper,optionsRm,options)
arguments
A BasisFunc;
Rc {Spin.mustBeSize(Rc,[3 3])}= diag([1 1 1]);
Rf {Spin.mustBeSize(Rf,[3 3])}= diag([1 1 1]);
tf {Spin.mustBeSize(tf,[3 1;1 3])}= ([0 0 0]);
optionsConvection.rightorleft {mustBeMember(optionsConvection.rightorleft,{'right','left'})}= 'right';
optionsOper.Oper = [];
optionsRm.Rm = POSCAR_read;
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
rightorleft = optionsConvection.rightorleft ;
optionsCell = namedargs2cell(options);
if ~isempty(optionsOper.Oper)
Rf = optionsOper.Oper.Rf;
if isempty(Rf)
optionsOper.Oper =optionsOper.Oper.attachRm(optionsRm.Rm);
end
Rf = optionsOper.Oper.Rf;
Rc = optionsOper.Oper.R;
tf = optionsOper.Oper.tf;
options.conjugate = optionsOper.Oper.conjugate;
options.antisymmetry = optionsOper.Oper.antisymmetry;
optionsCell = namedargs2cell(options);
Am = rotate(A,Rc,Rf,tf,rightorleft,optionsCell{:});
else
Am = rotate(A,Rc,Rf,tf,rightorleft,'Oper',optionsOper.Oper,optionsCell{:});
end
U  = InnerProduct(Am,A,'sym',options.sym);
end
