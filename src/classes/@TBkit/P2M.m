function M  = P2M(P)
            P_index = abs(P);
            P_phase = P./P_index;%exp(1i*angle(P));
            M  = zeros(length(P));
            for i =1:numel(P_index)
                M(i,P_index(i)) = P_phase(i);
            end
        end