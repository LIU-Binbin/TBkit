function SzM = Sz(SpinObj,options)
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
    SzM_L = repmat(S,[nSpin-1 1]);
else
    S = SpinObj.J;
    SzM_L = repmat(S,[nSpin-1 1]);
end
if options.full
    for a =1:nSpin
        SzM_L(a) = 2*(S+1-a)/2;
    end
    SzM = diag(SzM_L);
else
    SzM = InnerProduct(SpinObj,SzOper(SpinObj),'sym',options.sym );
end
end
