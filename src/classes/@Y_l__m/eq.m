function C = eq(A,B,options)
arguments
A Y_l__m;
B Y_l__m;
options.strict = false;
end
C = true;
if length(A) == 1 && length(B) == 1
if A.l ~= B.l
C = false;
return;
end
if A.m ~= B.m
C = false;
return;
end
if A.n ~= B.n
C = false;
return;
end
if options.strict
if A.coe ~= B.coe
C = false;
return;
end
end
elseif isvector(A) && isvector(B)
A = contract(A);
B = contract(B);
if length(A) ~= length(B)
C = false;return;
end
if isequal([A(:).l],[B(:).l])
C = false;return;
end
if isequal([A(:).m],[B(:).m])
C = false;return;
end
if isequal([A(:).n],[B(:).n])
C = false;return;
end
CoeLtmp = [A(:).coe;B(:).coe];
if rank(CoeLtmp) == 2
C = false;return;
end
if options.strict
if isequal(CoeLtmp(1,:),CoeLtmp(2,:))
C = false;
return;
end
end
else
C = false;
end
end
