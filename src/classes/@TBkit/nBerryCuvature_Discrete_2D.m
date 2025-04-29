function BC = nBerryCuvature_Discrete_2D(VV,Vk1,Vk2,Vk1k2)
            % Consider Vector them for speed up
            % VV    W_d_0;
            % Vk1   W_d_1;  略偏离kx的波函数
            % Vk2   W_d_2;  略偏离ky的波函数
            % Vk1k2 W_d_12; 略偏离kx，ky的波函数
            % ----- approach 1 -----
            % ----- approach 2 -----
            U = VV'*(Vk1*Vk1')*(Vk1k2*Vk1k2')*(Vk2*Vk2')*VV;
            BC = angle(eig(U));
        end