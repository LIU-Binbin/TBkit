function WAVEFuncL = smoothDegen(WAVEFuncL,V)
            % Degenerate perturbation theory
            [A,U] = eig(WAVEFuncL'*V*WAVEFuncL);
            %disp(WAVEFuncL)
            [A,~] = park.sorteig(U,A);
            WAVEFuncL = WAVEFuncL*A;
        end