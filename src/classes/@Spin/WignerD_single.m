function Matelement = WignerD_single(SpinObj1, SpinObj2, abc, rightorleft)
arguments
SpinObj1 Spin;
SpinObj2 Spin;
abc;
rightorleft = 'right';
end
assert(isrow(SpinObj1));
assert(isrow(SpinObj2));
alpha = sym(abc(1));
beta = sym(abc(2));
gamma = sym(abc(3));
Matelement = sym(0);
if strcmp(rightorleft, 'right')
for i = 1:length(SpinObj1)
for j = 1:length(SpinObj2)
if SpinObj1.J ~= SpinObj2.J
continue;
end
m1 = SpinObj1.Jz;
m2 = SpinObj2.Jz;
WignerD_single_element = sym(Y_l__m.d(SpinObj1.J, m1, m2, beta));
Matelement = Matelement + exp(1i * m1 * alpha) * WignerD_single_element * exp(1i * m2 * gamma);
end
end
end
end
