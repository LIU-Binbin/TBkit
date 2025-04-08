function dispAll(BasisFunction)
for i = 1:size(BasisFunction,1)
fprintf('============================================\n')
%fprintf('The %d / %d th BasisFunction:\n',i,numel(BasisFunction));
fprintf(string(BasisFunction(i,1).coe));
fprintf(string(BasisFunction(i,1).BFuncL));
fprintf(string(BasisFunction(i,1).spin));
for j = 2:size(BasisFunction,2)
fprintf(' + ');
fprintf(string(BasisFunction(i,j).coe));
fprintf(string(BasisFunction(i,j).BFuncL));
fprintf(string(BasisFunction(i,j).spin));
end
fprintf('\n')
end
fprintf('============================================\n')
end
