function A = uminus(A)
for i = 1:numel(A)
A(i).coe = -A(i).coe;
end
end
