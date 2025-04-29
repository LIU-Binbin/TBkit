function symPage = subspage(symPage,lhs,rhs)
            for i = 1:numel(lhs)
                symPage = subs(symPage,lhs(i),rhs(i));
            end
        end