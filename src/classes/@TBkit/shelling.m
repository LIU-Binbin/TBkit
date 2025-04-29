function SymVar = shelling(SymVar)
            % cos sin exp
            StrVar = char(SymVar);
            StrVar(end) = [];
            StrVar(1:4) = [];
            SymVar = str2sym(StrVar);
        end