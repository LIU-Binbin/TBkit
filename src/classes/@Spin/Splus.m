function Splus = Splus(SpinObj,options)
arguments
SpinObj Spin;
options.full = true;
options.sym = true;
end
optionsCell = namedargs2cell(options);
if options.full
Splus = Sx(SpinObj,optionsCell{:}) + 1i*Sy(SpinObj,optionsCell{:});
else
Splus = InnerProduct(SpinObj, SplusOper(SpinObj), 'sym', options.sym );
end
end
