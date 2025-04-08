function Lminus = Lminus(SpinObj,options)
arguments
SpinObj Spin;
options.full = true;
options.sym = true;
end
optionsCell = namedargs2cell(options);
Lminus = Lx(SpinObj,optionsCell{:}) - 1i*Ly(SpinObj,optionsCell{:});
end
