function OutExpr = nlm2atomic(l,m,n,options)
arguments
l = 1;
m = 0;
n = 0;
options.outputformat = 'sym';
end
L_dist = containers.Map([0,1,2,3,4,5,6],["s","p","d","f","g","h","i"]);
if strcmp(options.outputformat ,'sym')
if n ==0
str1 = "";
str2 = "";
else
str1 = num2str(n);
str2 = num2str(n);
end
if l >= 0 && l <=6
str1 = str1 + L_dist(l)+"_" + Y_l__m.lm2str(l,abs(m));
str2 = str2 + L_dist(l)+"_" + Y_l__m.lm2str(l,-abs(m));
elseif l < 0
else
end
if m == 0
OutExpr = sym(str1,'real');
else
if m < 0
OutExpr = 1/sqrt(2) * (sym(str1,'real')-1i*sym(str2,'real'));
else
OutExpr = (-1)^m/sqrt(2) * (sym(str1,'real')+1i*sym(str2,'real'));
end
end
end
end
