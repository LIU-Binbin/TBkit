function angle = name_angle(theta, Latex)
            arguments
                theta double ;
                Latex logical = false;
            end
            frac = rat(theta / pi, 1e-4);
            frac = simplify(str2sym(frac))*pi;
            if Latex
                angle = latex(frac);
            else
                angle = string(frac);
            end
        end