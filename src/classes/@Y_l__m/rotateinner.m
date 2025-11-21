function Aj = rotateinner(A,abc,RightorLeft,immproper,conjugate,antisymmetry)
arguments
    A
    abc
    RightorLeft
    immproper = false;
    conjugate = false;
    antisymmetry = false;
end

alpha = abc(1);
beta  = abc(2);
gamma = abc(3);

% ---------------------------
% 初始化输出 (保持你原逻辑)
% ---------------------------
Aj = A(1).HollowMe;

% ============================================================
% 循环对 A 中每一个 Y_l__m(ℓ,m) 进行复球谐旋转 (原逻辑)
% ============================================================
for j = 1:numel(A)

    A_j = A(j);         % 当前复球谐
    l   = A_j.l;
    m1  = A_j.m;
    c0  = A_j.coe;

    % 所有 Y_l__m(l,m2)
    A_L = Y_l__m(l);
    Ak  = A_j.HollowMe;

    % --------------------------------------------------------
    % 遍历所有 m2，对 Y_l__m(l,m1) → Y_l__m(l,m2) 构造 D^(l)
    % --------------------------------------------------------
    for i = 1:length(A_L)

        Ai = A_L(i);    % Y_l__m(l,m2)
        m2 = Ai.m;

        if l == 0
            coeff = 1;
        else
            % ----------------------------------------------
            % *** 核心修复：Wigner d 单元素必须严格用 m1,m2 顺序 ***
            %     并处理 beta = ±π 等奇点
            % ----------------------------------------------
            d_l = Y_l__m.d(l, m1, m2, beta);

            % 避免 ±pi 下出现 -0，+0，NaN，符号错乱
            if ~isfinite(d_l)
                d_l = 0;
            elseif abs(d_l) < 1e-14
                d_l = 0;
            end

            % Euler 角 (右/左旋) 一致
            coeff = exp(1i*RightorLeft*m1*alpha) * ...
                   d_l * ...
                   exp(1i*RightorLeft*m2*gamma);

            % 代码规定：复球谐旋转必须共轭 适应Convention
            coeff = conj(coeff);
        end
        coeff = c0*coeff;
        % improper 处理（你的原逻辑）
        if immproper
            coeff = (-1)^l * coeff;
        end

        % 外层 conjugate 选项（保持完全一致）
        if ~conjugate
            coeff = coeff;
        else
            coeff = conj( coeff);
        end


        Ai.coe = coeff;
        Ak = [Ak, Ai];
    end

    % 累积
    Aj = [Aj, Ak];
end

% 合并重复行（保持你原逻辑）
Aj = Aj.contractrow();

end


% function Aj = rotateinner(A,abc,RightorLeft,immproper,conjugate,antisymmetry)
% arguments
%     A
%     abc
%     RightorLeft
%     immproper = false;
%     conjugate = false;
%     antisymmetry = false;
% end
% alpha = (abc(1));
% beta = (abc(2));
% gamma = (abc(3));
% Aj = A(1).HollowMe;
% for j = 1:numel(A)
%     A_j = A(j);
%     A_L = Y_l__m(A_j.l);
%     Ak = A_j.HollowMe;
%     for i =  1:length(A_L)
%         Ai = A_L(i);
%         if A_j.l == 0
%             Ai.coe = 1;
%         else
%             m1 = A_j.m;
%             m2 = Ai.m;
%             WignerD_single_element = (Y_l__m.d(A_j.l,m1,m2,beta));
%             Ai.coe = Ai.coe*exp(1i*RightorLeft*m1*alpha)*WignerD_single_element*exp(1i*RightorLeft*m2*gamma);
%             Ai.coe = conj(Ai.coe);
%         end
%         if immproper
%             Ai.coe = (-1)^(A_j.l)*Ai.coe;
%         end
%         if ~conjugate
%             Ai.coe = Ai.coe * A_j.coe;
%         else
%             Ai.coe = conj(Ai.coe * A_j.coe);
%         end
%         Ak = [Ak,Ai];
%     end
%     Aj = [Aj.Ak];
% end
% Aj = Aj.contractrow();
% end
