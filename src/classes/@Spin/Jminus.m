function Jminus = Jminus(SpinObj,options)
arguments
SpinObj Spin;
options.full = true;
options.sym = true;
end
optionsCell = namedargs2cell(options);
Jminus = JX(SpinObj,optionsCell{:}) - 1i*JY(SpinObj,optionsCell{:});
end
