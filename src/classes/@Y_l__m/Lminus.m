function LminusMat = Lminus(Y_l__mObj,options)
arguments
Y_l__mObj Y_l__m;
options.full = false;
options.sym = true;
end
% optionsCell = namedargs2cell(options);

LminusMat = InnerProduct(Y_l__mObj, LminusOper(Y_l__mObj), 'sym', options.sym );

end