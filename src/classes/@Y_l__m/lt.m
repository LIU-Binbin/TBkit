function C = lt(A,B,options)
arguments
A Y_l__m;
B Y_l__m;
options.strict = false;
end
C = true;
if length(A) == 1 && length(B) == 1
if A.n > B.n
C = false;return;
end
if A.l > B.l
C = false;return;
end
if A.m > B.m
C = false;return;
end
if options.strict
if A == B
if A.coe >= B.coe
C = false;return;
end
end
end
else
C = length(A) < length(B);
end
end
