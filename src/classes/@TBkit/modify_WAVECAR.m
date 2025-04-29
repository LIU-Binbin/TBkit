function [WAVECAR_loop_] = modify_WAVECAR(WAVECAR_loop,BF_WAVECAR)
            WAVECAR_loop_ = zeros(size(WAVECAR_loop,1),size(BF_WAVECAR,2),size(BF_WAVECAR,3));
            for i = 1:size(WAVECAR_loop,3) % For k evolution 
                uk = WAVECAR_loop(:,:,i);
                vk = BF_WAVECAR(:,:,i);
                for j= 1:size(BF_WAVECAR,2) % for all occupied WAN band
                    WAVECAR_loop_(:,j,i) = sum(uk.*vk(:,j).',2);
                end
            end
        end