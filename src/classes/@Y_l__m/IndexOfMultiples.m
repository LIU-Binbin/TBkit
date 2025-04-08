function Ind = IndexOfMultiples(A)
T   = true(size(A));
off = false;
A   = A(:);
for iA = 1:numel(A)
if T(iA)
d = (A(iA) == A);
if sum(d) > 1
T(d) = off;
end
end
end
Ind = find(~T);
end
