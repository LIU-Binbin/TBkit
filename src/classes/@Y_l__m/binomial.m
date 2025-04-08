function c = binomial(n, k)
error(nargchk(1, 2, nargin));
if nargin == 1
validateattributes(n, {'numeric'}, {'scalar', 'real', 'nonnegative'});
c = diag(fliplr(pascal(floor(n) + 1))).';
else
assert(isscalar(n) || isscalar(k) || isequal(size(n),size(k)), ...
'Non-scalar arguments must have the same size.')
validateattributes([n(:); k(:)], {'numeric'}, {'real', 'nonnegative'});
c = exp(gammaln(n+1) - gammaln(k+1) - gammaln(n-k+1));
i = ( n==floor(n+.5) & k==floor(k+.5) );
c(i) = floor(c(i)+.5);
end
end
