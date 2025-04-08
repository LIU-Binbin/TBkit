function LzM = Lz(SpinObj,options)
arguments
SpinObj Spin;
options.full = true;
options.sym = true;
end
optionsCell = namedargs2cell(options);
LzM = Sz(SpinObj,optionsCell{:});
end
