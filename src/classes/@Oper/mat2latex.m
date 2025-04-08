function str = mat2latex(mat,accuracy)
mat = roundn(mat,-accuracy-1);
Size = size(mat);
str =  "\left(\begin{array}{";
for j =1:Size(2)
str = str+ 'c';
end
str = str+ '} ';
for i = 1:Size(1)
str = str + Oper.num2latex(mat(i,1),accuracy);
for j = 2:Size(2)
str = str +'&'+Oper.num2latex(mat(i,j),accuracy);
end
str = str+ '\\';
end
str = str + '\end{array}\right)';
end
