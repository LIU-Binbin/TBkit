function Y_l__mObj = LplusOper(Y_l__mobj)
Y_l__mObj = Y_l__mobj;
for i = 1:numel(Y_l__mobj)
    m = Y_l__mobj(i).m;
    l = Y_l__mobj(i).l;
    Y_l__mObj(i).coe = Y_l__mobj(i).coe * sqrt(l*(l+1)-m*(m+1));
    if m +1 <= l
        Y_l__mObj(i).m = Y_l__mobj(i).m + 1;
    else
        Y_l__mObj(i).m = Y_l__mobj(i).m ;
    end
end
% Y_l__mObj = Y_l__mObj.contract();
end