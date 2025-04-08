function disp(SpinObj)
strlist = string(SpinObj);
for i = 1:size(strlist,1)
fprintf(strlist(i,1));
end
end
