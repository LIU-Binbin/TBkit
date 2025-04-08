function phi = braket(phi1,Oper,phi2)
switch Oper
case 'S+'
S = phi2(1);
m = phi2(2);
phi = sqrt((S-m)*(S+m+1));
if m + 1 <= S
phi2(2) = m+1;
else
phi2(2) = -S;
end
if isequal(phi1,phi2)
else
phi= 0;
end
case 'S-'
S = phi2(1);
m = phi2(2);
phi = sqrt((S+m)*(S-m+1));
if m - 1 >= -S
phi2(2) = m-1;
else
phi2(2) = S;
end
if isequal(phi1,phi2)
else
phi= 0;
end
case '^S+'
S = phi2(1);
m = phi2(2);
phi = sqrt((S-m)*(S+m+1));
if m + 1 <= S
phi2(2) = S;
else
phi2(2) = -S;
end
if isequal(phi1,phi2)
else
phi= 0;
end
case '^S-'
S = phi2(1);
m = phi2(2);
phi = sqrt((S+m)*(S-m+1));
if m - 1 >= -S
phi2(2) = -S;
else
phi2(2) = S;
end
if isequal(phi1,phi2)
else
phi= 0;
end
end
end
