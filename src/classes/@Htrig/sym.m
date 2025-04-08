function Htrig_sym = sym(H_htrig,options)
arguments
H_htrig Htrig;
options.simple = false;
options.Type {mustBeMember(options.Type,{'exp','sincos','sin','cos','cot','tan','sinh','cosh','cosh','tanh',''})}= '';
end
Htrig_sym = H_htrig.Htrig_sym;
if strcmp(options.Type,'')
switch H_htrig.Type
case {'sincos'}
type = 'sincos';
case {'exp','list','mat'}
type = 'exp';
otherwise
type = 'sincos';
end
else
type = options.Type;
end
if options.simple
H_htrig = H_htrig.timtj_gen('sym');
Htrig_sym = simplify(Htrig_sym./H_htrig.timtj{3}.');
end
Htrig_sym = simplify(rewrite(Htrig_sym,type),'IgnoreAnalyticConstraints',true);
end
