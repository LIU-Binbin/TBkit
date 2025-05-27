function LzMat = Lz(Y_l__mObj,options)
arguments
Y_l__mObj Y_l__m;
options.full = false;
options.sym = true;
end
% optionsCell = namedargs2cell(options);

LzMat = InnerProduct(Y_l__mObj, LzOper(Y_l__mObj), 'sym', options.sym );

end