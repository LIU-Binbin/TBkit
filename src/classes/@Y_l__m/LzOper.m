function Y_l__mObj = LzOper(Y_l__mobj)
Y_l__mObj = Y_l__mobj;
for i = 1:numel(Y_l__mobj)
    m = Y_l__mobj(i).m;
    Y_l__mObj(i).coe = Y_l__mobj(i).coe * m;
end
% Y_l__mObj = Y_l__mObj;
end