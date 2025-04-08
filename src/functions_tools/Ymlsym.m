function orb_sym = Ymlsym(l, m, orb_symbk)
    % Function to generate orbital symmetry expressions based on l (angular quantum number) and m (magnetic quantum number)
    % orb_symbk is a backup if both l and m are complex (1i)
    %
    % Inputs:
    %   l - Angular quantum number (e.g., 0, 1, 2, 3, etc.)
    %   m - Magnetic quantum number (e.g., -l, ..., l)
    %   orb_symbk - Backup orbital symmetry if l and m are complex
    % Outputs:
    %   orb_sym - Symbolic representation of orbital symmetry

    % Return backup symmetry if l and m are imaginary
    if l == 1i && m == 1i
        orb_sym = orb_symbk;
        return;
    end

    % Switch case for different values of l
    switch l
        case 0
            % s-orbital (l=0), m=0
            if m == 0
                orb_sym = str2sym('1');
            else
                warning('Invalid input: check your POSCAR for l=0, m=%d', m);
            end
            
        case 1
            % p-orbital (l=1)
            switch m
                case 0
                    orb_sym = str2sym('z');
                case -1
                    orb_sym = str2sym('y');
                case 1
                    orb_sym = str2sym('x');
                otherwise
                    warning('Invalid input: check your POSCAR for l=1, m=%d', m);
            end
            
        case 2
            % d-orbital (l=2)
            switch m
                case 0
                    orb_sym = str2sym('z^2');
                case -1
                    orb_sym = str2sym('y*z');
                case 1
                    orb_sym = str2sym('x*z');
                case -2
                    orb_sym = str2sym('x*y');
                case 2
                    orb_sym = str2sym('x^2-y^2');
                otherwise
                    warning('Invalid input: check your POSCAR for l=2, m=%d', m);
            end

        case 3
            % f-orbital (l=3)
            switch m
                case 0
                    orb_sym = str2sym('z^3');
                case -1
                    orb_sym = str2sym('y*z^2');
                case 1
                    orb_sym = str2sym('x*z^2');
                case -2
                    orb_sym = str2sym('x*y*z');
                case 2
                    orb_sym = str2sym('(x^2-y^2)*z');
                case -3
                    orb_sym = str2sym('y*(3*x^2-y^2)');
                case 3
                    orb_sym = str2sym('x*(x^2-3*y^2)');
                otherwise
                    warning('Invalid input: check your POSCAR for l=3, m=%d', m);
            end

        % Higher orbital quantum numbers (l = -1, -2, -3, etc.)
        case {-1, -2, -3, -4, -5}
            orb_sym = handleHigherOrbitalSymmetries(l, m);
            
        otherwise
            warning('Invalid input: check your POSCAR for l=%d', l);
    end
end

% Helper function to handle higher orbital symmetries
function orb_sym = handleHigherOrbitalSymmetries(l, m)
    switch l
        case -1
            switch m
                case 1
                    orb_sym = str2sym('1+x');
                case 2
                    orb_sym = str2sym('1-x');
                otherwise
                    warning('Invalid input: check your POSCAR for l=-1, m=%d', m);
            end

        case -2
            switch m
                case 1
                    orb_sym = str2sym('3^(-1/2)-6^(-1/2)*x+2^(-1/2)*y');
                case 2
                    orb_sym = str2sym('3^(-1/2)-6^(-1/2)*x-2^(-1/2)*y');
                case 3
                    orb_sym = str2sym('3^(-1/2)+6^(-1/2)*x+6^(-1/2)*x');
                otherwise
                    warning('Invalid input: check your POSCAR for l=-2, m=%d', m);
            end

        case -3
            switch m
                case 1
                    orb_sym = str2sym('1+x+y+z');
                case 2
                    orb_sym = str2sym('1+x-y-z');
                case 3
                    orb_sym = str2sym('1-x+y-z');
                case 4
                    orb_sym = str2sym('1-x-y+z');
                otherwise
                    warning('Invalid input: check your POSCAR for l=-3, m=%d', m);
            end

        case -4
            switch m
                case 1
                    orb_sym = str2sym('3^(-1/2)-6^(-1/2)*x+2^(-1/2)*y');
                case 2
                    orb_sym = str2sym('3^(-1/2)-6^(-1/2)*x-2^(-1/2)*y');
                case 3
                    orb_sym = str2sym('3^(-1/2)+6^(-1/2)*x+6^(-1/2)*x');
                case 4
                    orb_sym = str2sym('2^(-1/2)*z+2^(-1/2)*z^2');
                case 5
                    orb_sym = str2sym('-2^(-1/2)*z+2^(-1/2)*z^2');
                otherwise
                    warning('Invalid input: check your POSCAR for l=-4, m=%d', m);
            end

        case -5
            switch m
                case 1
                    orb_sym = str2sym('6^(-1/2)-2^(-1/2)*x-12^(-1/2)*z^2+2^(-1)*(x^2-y^2)');
                case 2
                    orb_sym = str2sym('6^(-1/2)+2^(-1/2)*x-12^(-1/2)*z^2+2^(-1)*(x^2-y^2)');
                case 3
                    orb_sym = str2sym('6^(-1/2)-2^(-1/2)*x-12^(-1/2)*z^2-2^(-1)*(x^2-y^2)');
                case 4
                    orb_sym = str2sym('6^(-1/2)+2^(-1/2)*x-12^(-1/2)*z^2-2^(-1)*(x^2-y^2)');
                case 5
                    orb_sym = str2sym('6^(-1/2)-2^(-1/2)*z+3^(-1)*(z^2)');
                case 6
                    orb_sym = str2sym('6^(-1/2)+2^(-1/2)*z+3^(-1)*(z^2)');
                otherwise
                    warning('Invalid input: check your POSCAR for l=-5, m=%d', m);
            end

        otherwise
            warning('Invalid input: check your POSCAR for l=%d', l);
    end
end
