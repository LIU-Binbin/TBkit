function mustBeHalfInteger(A)
if ~isnumeric(A) && ~islogical(A)
throwAsCaller(createValidatorException('MATLAB:validators:mustBeNumericOrLogical'));
end
if ~isreal(A)
throwAsCaller(createValidatorException('MATLAB:validators:mustBeReal'));
end
A = 2*A;
if ~all(isfinite(A), 'all') || ~all(A == floor(A), 'all')
throwAsCaller(createValidatorException('MATLAB:validators:mustBeHalfInteger'));
end
end
