function name= pretty(SymOper,options)
arguments
    SymOper Oper ;
    options.full = true;
    options.Latex = false;
end
diminsion = size(SymOper.R,1);
switch diminsion
    case 1
        if SymOper.R(1,1) == 1
            rot_name = '1';
        else
            rot_name = 'I';
        end
    case 2
        if Oper.isclose(det(SymOper.R),1)
            theta = atan2(SymOper.R(1,2),SymOper.R(1,1));
            if Oper.isclose(theta,0)
                rot_name = '1';
            else
                if options.Latex
                    rot_name = "R\left("+Oper.name_angle(theta,options.Latex)+"\right)";
                else
                    rot_name = "R("+Oper.name_angle(theta)+")";
                end
            end
        elseif Oper.isclose(det(SymOper.R),-1)
            [val, vec] = eig(SymOper.R);
            [val,vec] =  Oper.sorteig(vec,val);
            if SymOper.allclose(diag(vec).',[-1,1])
                n = val(:,1).';
                if options.Latex
                    rot_name = "M\left("+mat2str(Oper.round_axis(n))+"\right)";
                else
                    rot_name = "M("+""+mat2str(Oper.round_axis(n))+")";
                end
            else
                error('?');
            end
        else
            error('det ~= \pm 1!');
        end
    case 3
        if Oper.isclose((det((SymOper.R))),1)
            [n,theta] = Oper.Rotation2nTheta(SymOper.R);
            if Oper.isclose(theta,0)
                rot_name = '1';
            else
                if options.Latex
                    rot_name = "R\left("+Oper.name_angle(theta,options.Latex)+","+mat2str(Oper.round_axis(n))+"\right)";
                else
                    rot_name = "R("+Oper.name_angle(theta)+","+mat2str(Oper.round_axis(n))+")";
                end
            end
        elseif Oper.isclose((det((SymOper.R))),-1)
            [n,theta] = Oper.Rotation2nTheta(-SymOper.R);
            if Oper.isclose(theta,0)
                rot_name = 'I';
            elseif Oper.isclose(theta,pi)
                if options.Latex
                    rot_name = "M\left("+mat2str(Oper.round_axis(n))+"\right)";
                else
                    rot_name = "M("+""+mat2str(Oper.round_axis(n))+")";
                end
            else
                if options.Latex
                    rot_name = "S\left("+Oper.name_angle(theta,latex)+","+mat2str(Oper.round_axis(n))+"\right)";
                else
                    rot_name = "S("+Oper.name_angle(theta)+","+mat2str(Oper.round_axis(n))+")";
                end
            end
        else
            error('det ~= 1!');
        end
end
if options.full
    if options.Latex
        name = "\begin{aligned} U H(\mathbf{{k}})";
        if SymOper.conjugate
            name = name + "^*";
        end
        name = name + " U^{{-1}} &= ";
        if SymOper.antisymmetry
            name = name + "-";
        end
        name = name + "H(";
        if SymOper.conjugate
            name = name + "-";
        end
        name = name + "R\mathbf{{k}}) \\";
    else
        name = "U·H(k)";
        if SymOper.conjugate
            name = name + "*";
        end
        name = name + "·U^-1 = ";
        if SymOper.antisymmetry
            name = name + "-";
        end
        name = name + "H(";
        if SymOper.conjugate
            name = name + "-";
        end
        name = name + "R·k)\n";
    end
    if options.Latex
        name = name + 'R &= '+rot_name + '\\';
    else
        name = name + 'R = '+rot_name + '\n';
    end
    if ~isnan(SymOper.U)
        if isa(SymOper.U,'sym')
            if options.Latex
                Umat = latex(SymOper.U);
                name = name + "U &= "+Umat+'\end{aligned}';
            else
                name = name + "U = "+mat2str(string(SymOper.U))+'\n';
            end
        else
            if options.Latex
                Umat = Oper.mat2latex(SymOper.U,6);
                name = name + "U &= "+Umat+'\end{aligned}';
            else
                name = name + "U = "+mat2str(roundn(SymOper.U,-4),3)+'\n';
            end
        end
    else
        if options.Latex
            name = name + '\end{aligned}';
        else
            name = name + '\n';
        end
    end
else
    if SymOper.conjugate && ~SymOper.antisymmetry
        if options.Latex
            az_name = "\mathcal{T}";
        else
            az_name = "T";
        end
    elseif SymOper.conjugate && SymOper.antisymmetry
        if options.Latex
            az_name = "\mathcal{P}";
        else
            az_name = "P";
        end
    elseif ~SymOper.conjugate && SymOper.antisymmetry
        if options.Latex
            az_name = "\mathcal{C}";
        else
            az_name = "C";
        end
    else
        az_name = "";
    end
    if strcmp(rot_name ,'1') && ~strcmp(az_name,"")
        if ~strcmp(az_name,"")
            name = az_name + " "+az_name;
        else
            name = az_name + ""+az_name;
        end
    else
        if ~strcmp(az_name,"")
            name = rot_name + " "+az_name;
        else
            name = rot_name + ""+az_name;
        end
    end
end
end
