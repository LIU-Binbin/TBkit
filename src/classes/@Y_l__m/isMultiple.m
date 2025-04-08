function T = isMultiple(A)
T        = false(size(A));
[S, idx] = sort(A(:).');
m        = [false, diff(S) == 0];
if any(m)
m(strfind(m, [false, true])) = true;
T(idx) = m;
end
end
