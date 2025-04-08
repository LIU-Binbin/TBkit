function group = square(tr, ph, generators, spin,options,propArgs)
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
if ~all(all(isnan(ph))) && ~all(all(all(isnan(spin))))
raise ValueError('If `ph` is True, `spin` may not be provided, as it is not '...
'possible to deduce the unitary representation of particle-hole symmetry '...
'from spin alone. In this case construct the particle-hole operator manually.');
end
Mx = Oper.mirror([1, 0],nan,spin,optionsCell{:},propertyCell{:});
C4 = Oper.rotation(1/4, nan,false,nan, spin,optionsCell{:},propertyCell{:});
gens = [Mx, C4];
if tr
propArgs.antisymmetry = true;
propertyCell = namedargs2cell(propArgs);
TR = Oper.time_reversal(2, spin,propertyCell{:});
gens = [gens,TR];
end
if ph
propArgs.antisymmetry = true;
propArgs.conjugate = true;
PH = Oper.particle_hole(2,propertyCell{:});
gens = [gens,PH];
end
if generators
group =  gens;
else
%fprintf('group muplicity: %d\n',gens.order);
group =  gens.generate_group();
end
end
