function SyM = Sy(SpinObj,options)
arguments
SpinObj Spin;
options.full = false;
options.sym = true;
end
if options.full
SpinObj = SpinObj.BasisGen();
end
optionsCell = namedargs2cell(options);
nSpin = length(SpinObj);
if options.full
if nSpin ~= SpinObj(1).J*2+1
error('!');
end
end
if options.sym
S = sym(SpinObj(1).J);
SyM_L = repmat(S,[nSpin-1 1]);
else
S = SpinObj.J;
SyM_L = repmat(S,[nSpin-1 1]);
end
if options.full
for a =1:nSpin-1
b = a+1;
SyM_L(a) = -1i*sqrt((S+1)*(a+b-1)-a*b)/2;
end
SyM = diag(SyM_L,1) + diag(conj(SyM_L),-1);
else
SyM = 1/2i*(Splus(SpinObj,optionsCell{:}) - Sminus(SpinObj,optionsCell{:}));
end
end
