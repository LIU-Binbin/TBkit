function JyM = JY(SpinObj,options)
arguments
SpinObj Spin;
options.full = true;
options.sym = true;
end
optionsCell = namedargs2cell(options);
JyM = Sy(SpinObj,optionsCell{:});
end
