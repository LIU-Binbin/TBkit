function B = contractrow(A,options)
arguments
A Y_l__m;
options.forgetcoe = false;
end
B = A;
A = cleanrow(A);
lmL = ([A(:).l;A(:).m;]).';
coeL = [A(:).coe];
[lmL_unique,~,~] = unique(lmL,'rows');
if isempty(lmL_unique)
return;
end
if options.forgetcoe
A = repmat(A(1),[1 length(lmL_unique)]);
for i =1:length(lmL_unique)
A(i).l = lmL_unique(i,1);
A(i).m = lmL_unique(i,2);
A(i).coe = 1;
end
else
B = repmat(A(1),[1 size(lmL_unique,1)]);
for i = 1:size(lmL_unique,1)
lml_tmp = lmL_unique(i,:);
[~,index_ic] = ismember(lmL,lml_tmp,'rows');
B(1,i).l = lmL_unique(i,1);
B(1,i).m = lmL_unique(i,2);
B(1,i).coe = sum(coeL(logical(index_ic)));
end
end
B = cleanrow(B);
end
