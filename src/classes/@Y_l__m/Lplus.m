function LplusMat = Lplus(Y_l__mObj,options)
arguments
Y_l__mObj Y_l__m;
options.full = false;
options.sym = true;
end
% optionsCell = namedargs2cell(options);

LplusMat = InnerProduct(Y_l__mObj, LplusOper(Y_l__mObj), 'sym', options.sym );
% bug
% LplusMat
end
