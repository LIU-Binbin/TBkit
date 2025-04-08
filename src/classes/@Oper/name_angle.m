function angle = name_angle(theta, Latex)
arguments
theta  ;
Latex logical = false;
end
if isa(theta,'sym')
frac = simplify(theta,'IgnoreAnalyticConstraints',true);
else
frac = rat(theta / pi, 1e-4);
frac = simplify(str2sym(frac))*pi;
end
if Latex
angle = latex(frac);
else
angle = string(frac);
end
end
