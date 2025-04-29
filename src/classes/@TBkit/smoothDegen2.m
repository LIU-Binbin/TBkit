function WAVEFuncL = smoothDegen2(WAVEFuncL)
            % Wave1 Todo
            WAVEFuncLabs =abs(WAVEFuncL(:,1,:));
            [~,maxWAVEFuncLabs] = max(WAVEFuncLabs,[],1);
            [~,SelectOrb] = max(sum(maxWAVEFuncLabs,2));
            for i = 1:size(WAVEFuncL,2)
                for k = 1:size(WAVEFuncL,3)
                    WAVEFuncLtmp = exp(-1i*angle(WAVEFuncL(SelectOrb,1,i)))*WAVEFuncL(:,:,i);
                    WAVEFuncL(:,:,i) = WAVEFuncLtmp;
                end
            end
            %WAVEFuncLimag =imag(WAVEFuncL);
        end