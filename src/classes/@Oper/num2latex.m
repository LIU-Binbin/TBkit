function str = num2latex(num,accuracy)
arguments
num double
accuracy double  =6;
end
digits(accuracy);
theta = angle(num);
r = abs(num);
theta_name = Oper.name_angle(theta,1);
theta_expr = "e^{i "+theta_name+"}";
if Oper.isclose(imag(num),0)
str = latex(sym(num));
elseif Oper.isclose(real(num),0)
str = latex(sym(num));
else
if Oper.isclose(r,1)
str = theta_expr;
else
str = latex(sym(r))+theta_expr;
end
end
end
