function labelcut_list = labelcut_list_gen(Atom_num)
            n = length(Atom_num);
            beginline =8 ;
            sum_n =beginline;
            for i =1:n
                sum_n2 = sum_n+Atom_num(i)-1;
                labelcut_list(i,:) = [sum_n sum_n2];
                sum_n = sum_n2+1;
            end

        end