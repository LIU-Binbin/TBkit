function orbsym_n = subs_xyz(orbsym,Rlmn)
            % subs ok
            %disp(Rlmn);
            if orbsym == sym(1)
                orbsym_n = 1;
            elseif orbsym == sym('x')
                orbsym_n = Rlmn(1);

            elseif orbsym == sym('y')
                orbsym_n = Rlmn(2);

            elseif orbsym == sym('z')
                orbsym_n = Rlmn(3);

            else
                orbsym_n = 1;
            end
            % l m n fixed ?
            % orbsym_n =abs(orbsym_n);
            %disp(orbsym_n);
        end