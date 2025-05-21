function SymOper_out = times(SymOper1,SymOper2)
m = length(SymOper1);
n = length(SymOper2);
if m == 1&& n ==1
    SymOper_out = SymOper1;
    if all(all((isnan(SymOper1.U))))
        SymOper_out.U = nan;
    elseif SymOper1.conjugate
        SymOper_out.U = SymOper1.U*conj(SymOper2.U);
    else
        SymOper_out.U = SymOper1.U*(SymOper2.U);
    end
    SymOper_out.R =  SymOper1.R*(SymOper2.R);
    SymOper_out.t =  SymOper2.t*SymOper1.R.'+SymOper1.t;
    SymOper_out.conjugate = xor(SymOper1.conjugate,SymOper2.conjugate);
    SymOper_out.antisymmetry = xor(SymOper1.antisymmetry,SymOper2.antisymmetry);
    SymOper_out.strict_eq = SymOper1.strict_eq || SymOper2.strict_eq;
else
    error('.* is elemenary operator');
end
end
