function NlmExpr = Nlm(l,m,seed,options)
arguments
    l {mustBeNonnegative,mustBeInteger};
    m {mustBeInteger}
    seed sym{} = (sym('theta','real'));
    options.ClosedForm = false;
    options.triangle = true;
end
optionsCell = namedargs2cell(options);
NlmExpr = (-1)^m * sqrt(((l+1/2))*factorial((l-m))/factorial((l+m)))* ...
    Y_l__m.Plm(l,m,seed,optionsCell{:});
end
