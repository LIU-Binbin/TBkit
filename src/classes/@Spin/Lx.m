function LxM = Lx(SpinObj,options)
arguments
SpinObj Spin;
options.full = true;
options.sym = true;
end
optionsCell = namedargs2cell(options);
LxM = Sx(SpinObj,optionsCell{:});
end
