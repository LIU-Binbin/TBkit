function expression = formula(YlmObj,options)
arguments
    YlmObj Y_l__m;
    options.vpa = true;
    options.explicit = false;
    options.cart = false;
end
expression =sym(0);
optionsCell = namedargs2cell(options);
if options.cart
    expression = repmat(expression,[size(YlmObj,1),1]);
    for j = 1:size(YlmObj,1)
        for i = 1:size(YlmObj,2)
            if ~YlmObj(j,i).hollow
                expression(j) = expression(j)  + YlmObj(j,i).coe*Y_l__m.nlm2atomic(YlmObj(j,i).l,YlmObj(j,i).m,YlmObj(j,i).n,'outputformat','sym');
            end
        end
    end
    expression = simplify(expression);
    return;
end
if options.explicit
    for i = 1:length(YlmObj)
        expression = expression +YlmObj(i).coe*explicitformula(YlmObj(i),optionsCell{:});
    end
else
    if options.vpa
        for i = 1:length(YlmObj)
            Symble = TBkit.SymbolicVarible('Y',YlmObj(i).m,YlmObj(i).l);
            expression = expression + vpa(YlmObj(i).coe)*Symble;
        end
    else
        for i = 1:length(YlmObj)
            Symble = TBkit.SymbolicVarible('Y',YlmObj(i).m,YlmObj(i).l);
            expression = expression + YlmObj(i).coe*Symble;
        end
    end
end
end
