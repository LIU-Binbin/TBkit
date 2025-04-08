function JzM = JZ(SpinObj,options)
arguments
SpinObj Spin;
options.full = true;
options.sym = true;
end
optionsCell = namedargs2cell(options);
JzM = Sz(SpinObj,optionsCell{:});
end
