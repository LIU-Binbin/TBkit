function n = round_axis(n)
if isa(n,'sym')
n = simplify(n,'IgnoreAnalyticConstraints',true);
n = string(n);
else
n = roundn(n, -2);
end
end
