function B = contractrow(A, options)
arguments
A Spin;
options.forgetcoe = false;
end
B = A;
A = cleanrow(A);
jmL = [A(:).J; A(:).Jz].';
coeL = [A(:).coe];
[jmL_unique, ~, ~] = unique(jmL, 'rows');
if isempty(jmL_unique)
return;
end
if options.forgetcoe
A = repmat(A(1), [1 length(jmL_unique)]);
for i = 1:length(jmL_unique)
A(i).J = jmL_unique(i, 1);
A(i).Jz = jmL_unique(i, 2);
A(i).coe = 1;
end
else
B = repmat(A(1), [1 size(jmL_unique, 1)]);
for i = 1:size(jmL_unique, 1)
jml_tmp = jmL_unique(i, :);
[~, index_ic] = ismember(jmL, jml_tmp, 'rows');
B(1, i).J = jmL_unique(i, 1);
B(1, i).Jz = jmL_unique(i, 2);
B(1, i).coe = sum(coeL(logical(index_ic)));
end
end
B = cleanrow(B);
end
