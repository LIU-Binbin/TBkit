function disp(YlmObj,options)
arguments
YlmObj Y_l__m;
options.vpa = true;
options.explicit = true;
options.cart = false
end
optionsCell = namedargs2cell(options);
if length(YlmObj)== 1
disp(YlmObj.explicitformula(optionsCell{:}));
else
disp(YlmObj.formula(optionsCell{:}));
end
end
