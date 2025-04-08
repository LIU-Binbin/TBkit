function SymOper = Ugen(SymOper,Basis,options)
arguments
SymOper Oper
Basis
options.Rm = POSCAR_read;
options.sym = false;
options.center = [0,0,0];
end
optionsCell = namedargs2cell(options);
for i =1:numel(SymOper)
if isempty(SymOper(i).Rf)
SymOper(i) = SymOper(i).attachRm(options.Rm);
end
if isnan(SymOper(i).U)
SymOper(i).U =Basis.rotation('Oper', SymOper(i),optionsCell{:});
end
end
end
