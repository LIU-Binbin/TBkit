function Kabc_ = Kabc(Aab, VEC_ki, dEnm, Delta_abc, eps)
%--------------------------------------------------------
% 计算三阶张量 K_{nm}^{abc}
%
% 输入:
%   Aab        - shift_vector 输出, 大小 (Nbands,Nbands,3,3)
%   VEC_ki     - 速度矩阵 v^a, 大小 (Nbands,Nbands,3)
%   dEnm       - 能量差 E_n - E_m, 大小 (Nbands,Nbands)
%   Delta_abc  - Delta_{nm}^a, 大小 (Nbands,Nbands,3)
%   eps        - 容差 (避免除零), 默认 1e-8
%
% 输出:
%   Kabc       - K_{nm}^{abc}, 大小 (Nbands,Nbands,3,3,3)
% Ref:
% $$ K_{nm}^{abc} = -\frac{v_{nm}^{a}}{\varepsilon_{nm}^{3}}\mathcal{A}_{mn;b}^{c}+i\frac{v_{nm}^{a}\Delta_{nm}^{b}v_{mn}^{c}}{\varepsilon_{nm}^{5}}, \quad \Delta_{nm}^a=v_n^a-v_m^a. $$
%--------------------------------------------------------

if nargin < 5
    eps = 1e-8;
end

Nbands = size(VEC_ki,1);
Kabc_ = zeros(Nbands,Nbands,3,3,3);

for a = 1:3
    for b = 1:3
        for c = 1:3
            for n = 1:Nbands
                for m = 1:Nbands
                    if n == m
                        continue;
                    end
                    Enm = dEnm(n,m);
                    if abs(Enm) < eps
                        continue;
                    end

                    % v_{nm}^a, v_{mn}^c
                    vnm_a = VEC_ki(n,m,a);
                    vmn_c = VEC_ki(m,n,c);

                    % Delta_{nm}^b
                    Delta_b = Delta_abc(n,m,b);

                    % A_{mn;b}^c  注意公式是 mn 索引
                    Amn_bc = Aab(m,n,b,c);

                    % 公式两部分
                    term1 = -(vnm_a / (Enm^3)) * Amn_bc;
                    term2 = 1i * (vnm_a * Delta_b * vmn_c) / (Enm^5);

                    Kabc_(n,m,a,b,c) = term1 + term2;
                end
            end
        end
    end
end

end
