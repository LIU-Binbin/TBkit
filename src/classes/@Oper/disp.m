function str = disp(SymOper,options)
arguments
SymOper Oper;
options.full logical =false;
end
len = length(SymOper);
if len == 1
builtin('disp',SymOper);
try
str = SymOper.pretty();
fprintf(str);
catch
end
return;
end
if options.full
for i =1:len
if len > 1
fprintf("=============== %dth Oper ===============\n",i);
end
builtin('disp',SymOper(i));
str = SymOper(i).pretty();
fprintf(str);
end
else
for i =1:len
str = SymOper(i).pretty('full',false);
fprintf("%3d : %s\n",i,str);
end
end
end
