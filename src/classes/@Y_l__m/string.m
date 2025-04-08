function str = string(YlmObj,options)
arguments
YlmObj Y_l__m;
options.vpa = true;
options.explicit = false;
options.cart = true;
end
optionsCell = namedargs2cell(options);
str = string(formula(YlmObj,optionsCell{:}));
end
