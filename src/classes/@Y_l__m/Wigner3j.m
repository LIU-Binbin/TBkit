function  W = Wigner3j( j123, m123 )
j1 = j123(1); j2 = j123(2); j3 = j123(3);
m1 = m123(1); m2 = m123(2); m3 = m123(3);
if any( j123 < 0 ),
error( 'The j must be non-negative' )
elseif any( rem( [j123, m123], 0.5 ) ),
error( 'All arguments must be integers or half-integers' )
elseif any( rem( (j123 - m123), 1 ) )
error( 'j123 and m123 do not match' );
end
if ( j3 > (j1 + j2) ) || ( j3 < abs(j1 - j2) ) ...
|| ( m1 + m2 + m3 ~= 0 ) ...
|| any( abs( m123 ) > j123 ),
W = 0;
return
end
if ~any( m123 ) && rem( sum( j123 ), 2 ),
W = 0;
return
end
t1 = j2 - m1 - j3;
t2 = j1 + m2 - j3;
t3 = j1 + j2 - j3;
t4 = j1 - m1;
t5 = j2 + m2;
tmin = max( 0,  max( t1, t2 ) );
tmax = min( t3, min( t4, t5 ) );
t = tmin : tmax;
W = sum( (-1).^t .* exp( -ones(1,6) * gammaln( [t; t-t1; t-t2; t3-t; t4-t; t5-t] +1 ) + ...
gammaln( [j1+j2+j3+1, j1+j2-j3, j1-j2+j3, -j1+j2+j3, j1+m1, j1-m1, j2+m2, j2-m2, j3+m3, j3-m3] +1 ) ...
* [-1; ones(9,1)] * 0.5 ) ) * (-1)^( j1-j2-m3 );
if isnan( W )
warning( 'MATLAB:Wigner3j:NaN', 'Wigner3J is NaN!' )
elseif isinf( W )
warning( 'MATLAB:Wigner3j:Inf', 'Wigner3J is Inf!' )
end
end
