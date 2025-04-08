function Outstr = l__m2str(l,m)
if l <1
end
switch l
case 0
switch m
case 0
Outstr = "s";
end
case 1
switch m
case 0
Outstr = "z";
case -1
Outstr = "y";
case 1
Outstr = "x";
end
case 2
switch m
case 0
Outstr = "z2";
case -1
Outstr = "yz";
case 1
Outstr = "xz";
case -2
Outstr = "xy";
case 2
Outstr = "x2my2";
end
case 3
switch m
case 0
Outstr = "z3";
case -1
Outstr = "yz2";
case 1
Outstr = "xz2";
case -2
Outstr = "xyz";
case 2
Outstr = "zx2my2";
case -3
Outstr = "3x2ymy3";
case 3
Outstr = "x2m3y2";
end
otherwise
Outstr = num2str(abs(m));
if sign(m) == -1
Outstr = [Outstr,'__bar'];
end
end
end
