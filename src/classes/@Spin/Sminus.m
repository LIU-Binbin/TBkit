function Sminus = Sminus(SpinObj,options)
arguments
SpinObj Spin;
options.full = true;
options.sym = true;
end
optionsCell = namedargs2cell(options);
if options.full
Sminus = Sx(SpinObj,optionsCell{:}) - 1i*Sy(SpinObj,optionsCell{:});
else
Sminus = InnerProduct(SpinObj, SminusOper(SpinObj), 'sym', options.sym );
end
end
