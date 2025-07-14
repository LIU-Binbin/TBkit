function zak_phase_tot = ZakPhase(HR, num_occ, kpath_file)
% ZakPhase 计算紧束缚模型在某一K路径方向上的Zak Phase（总 Berry Phase）
%
% 用法:
%   zak_phase_tot = ZakPhase(HR, num_occ, kpath_file)
%
% 输入参数:
%   HR         - HR类对象（如HR_up或HR_dn）
%   num_occ    - 占据带数量（整数）
%   kpath_file - K点路径文件名（如 'KPOINTS_BP'）
%
% 输出:
%   zak_phase_tot - 总的 Zak Phase，单位为弧度，范围 [0, 2π]
%
% 注意事项:
%   - 本函数保持轨道中心位置 HR.orbL 不变，仅使用相位因子进行 Berry phase 修正
%   - 所有占据带的 Berry 相位通过“Wilson loop”方法积累
%   - 默认使用 det(Uk' * Uk+1) 的积累方式，适合有能隙的绝缘体

    % 设置计算路径
    HR = HR < kpath_file;

    % 生成能带与波函数
    [EIGENCAR, WAVECAR] = HR.EIGENCAR_gen('LWAVE', true);

    % 提取维度信息
    [num_wan, ~, Nk] = size(WAVECAR);

    % 对每个k点的波函数进行cell-periodic波函数修正（exp(-i k·r)）
    for k = 1:Nk
        % 生成轨道位置对应的相位因子（列向量）
        phase_correction = exp(-1i * 2 * pi * (HR.orbL * HR.klist_frac(k,:)'));  % num_wan × 1
        % 应用于每一个占据态波函数
        for occ = 1:num_occ
            WAVECAR(:, occ, k) = WAVECAR(:, occ, k) .* phase_correction;
        end
    end

    % 初始化Zak phase总值
    theta = 0;

    % 主循环：计算 U_k^† · U_{k+1} 的行列式的角度，并累加
    for k = 1:Nk-1
        Uk   = WAVECAR(:,1:num_occ,k);     % num_wan × num_occ
        Ukp1 = WAVECAR(:,1:num_occ,k+1);   % num_wan × num_occ
        theta = theta + angle(det(Uk' * Ukp1));
    end

    % 闭合段：U_{Nk} 到 U_1
    Uk   = WAVECAR(:,1:num_occ,Nk);
    Ukp1 = WAVECAR(:,1:num_occ,1);
    theta = theta + angle(det(Uk' * Ukp1));

    % 最终 Zak Phase，规约到 [0, 2π]
    zak_phase_tot = mod(theta, 2*pi);
end
