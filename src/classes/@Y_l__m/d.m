function SymExpr = d(j,m1,m2,seed)
arguments
    j double{mustBePositive};
    m1 double
    m2 double
    seed sym
end
theta = seed;
if (m1) < (m2)
    SymExpr = (-1)^(m2-m1)*Y_l__m.d(j,m2,m1,seed);
    return;
end
if (m1 <=0 ||(m1>0 && abs(m2)> m1))&& m2 <0
    SymExpr = Y_l__m.d(j,-m2,-m1,seed);
    return;
end
if isinteger(j)
    if m1==0 && m2 == 0
        SymExpr = legendreP(j,cos(theta));
        return
    end
    if m1 < 0 && m2 == 0
        SymExpr = Y_l__m(j,m1,sqrt(4*sym(pi)/(2*j+1)),1).explicitformula(...
            'vpa',false,'seed',[char(seed),'0']);
        return;
    end
end
switch j
    case 1/2
        if m1 ==1/2 && m2 ==1/2
            SymExpr = cos(theta/2);
        elseif m1 ==1/2 && m2 ==-1/2
            SymExpr = -sin(theta/2);
        else
            SymExpr = sym(0);
        end
    case 1
        if m1 ==1 && m2 ==1
            SymExpr = 0.5*(1+cos(theta));
        elseif m1 ==1 && m2 ==-1
            SymExpr = 0.5*(1-cos(theta));
        elseif m1 ==1 && m2 ==0
            SymExpr = -1/2^0.5*sin(theta);
        elseif m1 ==0 && m2 ==0
            SymExpr = cos(theta);
        else
            SymExpr = sym(0);
        end
    case 3/2
        switch m1
            case 3/2
                switch m2
                    case 3/2
                        SymExpr = 1/2*(cos(theta)+1)*cos(theta/2);
                    case 1/2
                        SymExpr = -sqrt(3)/2*(cos(theta)+1)*sin(theta/2);
                    case -1/2
                        SymExpr = sqrt(3)/2*(-cos(theta)+1)*cos(theta/2);
                    case -3/2
                        SymExpr = -1/2*(-cos(theta)+1)*sin(theta/2);
                end
            case 1/2
                switch m2
                    case 1/2
                        SymExpr = 1/2*(3*cos(theta)-1)*cos(theta/2);
                    case -1/2
                        SymExpr = -1/2*(3*cos(theta)+1)*sin(theta/2);
                end
        end
    case 2
        switch m1
            case 2
                switch m2
                    case 2
                        SymExpr = 1/4*(cos(theta)+1)^2;
                    case 1
                        SymExpr = -1/2*sin(theta)*(cos(theta)+1);
                    case -1
                        SymExpr = sqrt(3)/2*(-cos(theta)+1)*cos(theta/2);
                    case -2
                        SymExpr = -1/2*sin(theta)*(-cos(theta)+1);
                    case 0
                        SymExpr = sqrt(3/8)*sin(theta)^2;
                end
            case 1
                switch m2
                    case 1
                        SymExpr = 1/2*(cos(2*theta)+cos(theta));
                    case -1
                        SymExpr = 1/2*(-cos(2*theta)+cos(theta));
                    case 0
                        SymExpr = -sqrt(3/8)*sin(2*theta);
                end
            case 0
                switch m2
                    case 0
                        SymExpr = 1/2*(3*cos(theta)^2-1);
                end
        end
    otherwise
        SymExpr = Y_l__m.d_mm__j(j,m1,m2,seed);
end
end
