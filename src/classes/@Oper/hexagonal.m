function group = hexagonal(tr, ph, generators, spin,options,propArgs)
arguments
tr logical =true;
ph logical =true;
generators logical = false;
spin double = nan;
options.sym = false;
propArgs.?Oper;
end
optionsCell = namedargs2cell(options);
propertyCell = namedargs2cell(propArgs);
if ~isnan(ph) && ~isnan(spin)
raise ValueError('If `ph` is True, `spin` may not be provided, as it is not '...
'possible to deduce the unitary representation of particle-hole symmetry '...
'from spin alone. In this case construct the particle-hole operator manually.');
end
if isnan(spin)
I = Oper.inversion(3,propertyCell{:});
else
U = Oper.spin_rotation(zeros(3), spin);
I = Oper.inversion(3, U,propertyCell{:});
end
C2x = Oper.rotation(1/2, [1, 0, 0],false,nan,spin,optionsCell{:},propertyCell{:});
C6 = Oper.rotation(1/6, [0, 0, 1],false,nan,spin,optionsCell{:},propertyCell{:});
gens = [I ,C2x,C6];
if tr
propArgs.antisymmetry = true;
propertyCell = namedargs2cell(propArgs);
TR = Oper.time_reversal(3, spin,propertyCell{:});
gens = [gens,TR];
end
if ph
propArgs.antisymmetry = true;
propArgs.conjugate = true;
propertyCell = namedargs2cell(propArgs);
PH = Oper.particle_hole(3,propertyCell{:});
gens = [gens,PH];
end
if generators
group =  gens;
else
%fprintf('group muplicity: %d\n',gens.order);
group =  gens.generate_group();
end
end
