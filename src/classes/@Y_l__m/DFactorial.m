function DFact = DFactorial(n)
if isempty(n) || (numel(n) > 1) || (n < 0) || (mod(n,1) ~= 0)
error('The sky is falling. n must be scalar, non-negative, integer.')
end
DFact = 1;
if n > 1
start = 1 + mod(n + 1,2);
DFact = prod(start:2:n);
end
end
