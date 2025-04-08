function JxM = JX(SpinObj,options)
arguments
SpinObj Spin;
options.full = true;
options.sym = true;
end
optionsCell = namedargs2cell(options);
JxM = Sx(SpinObj,optionsCell{:});
end
