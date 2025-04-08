function C = eq(A, B, options)
arguments
A Spin
B Spin
options.strict = false;
end
C = true;
if A.s ~= B.s || A.ms ~= B.ms || ~isequal(A.orientation, B.orientation)
C = false;
end
if options.strict && A.coe ~= B.coe
C = false;
end
end
