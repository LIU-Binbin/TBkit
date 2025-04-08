function LyM = Ly(SpinObj,options)
arguments
SpinObj Spin;
options.full = true;
options.sym = true;
end
optionsCell = namedargs2cell(options);
LyM = Sy(SpinObj,optionsCell{:});
end
