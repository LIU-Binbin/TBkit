function C = minus(A, B, options)
arguments
A
B
options.autocontract = true;
end
if isa(A, 'Spin') && isa(B, 'Spin')
C = A;
C = [C, -B];
if options.autocontract
C = contract(C);
end
end
end
