function l_num = orb2l(input)
input = string(input);
switch input
    case {'0','s','S'}
        l_num = 0;
    case {'1','p','P','px','py','pz','Px','Py','Pz','p_x','p_y','p_z'}
        l_num = 1;
    case {'2','d','D','dx2-y2','dz2','dxy','dyx','dxz','dyz','dzx','dzy','dx2y2'}
        l_num = 2;
    case {'3','f'}
        l_num = 3;
    case {'-1','sp'}
        l_num = -1;
    case {'-2','sp2'}
        l_num = -2;
    case {'-3','sp3'}
        l_num = -3;
    case {'-4','sp3d'}
        l_num = -4;
    case {'-5','sp3d2'}
        l_num = -5;
    otherwise
        l_num = 1i;
end


end