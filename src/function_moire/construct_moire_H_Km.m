
function Hm = construct_moire_H_Km(kpt, Hfun, Km, Gm, Vm, N)
% construct_moire_H_Km  在 Km 周围三个最近邻 BZ 构造 moiré 哈密顿量
%
%  输入：
%    kpt  - 1×3,   从 Km 出发的局部偏移（习惯上我们把 klist_r 视为相对 Km 的偏移）
%    Hfun - 函数句柄 @(kx,ky,kz) → 原始晶格的 N×N 哈密顿量
%    Km   - 1×3,   moiré K 点坐标
%    Gm   - 3×3,   三个 moiré 基矢（行向量）
%    Vm   - N×N,   Moiré 势（中心与任一邻点的耦合矩阵）
%    N    - 标量,  单个 block 的维度
%
%  输出：
%    Hm   - 4N×4N, 构造好的 moiré 哈密顿量，恢复 C3 对称

    % 预分配：先构造对角块
    % Block 1: H0 = Hfun(Km + kpt)
    [e0, psi0] = vasplib.EIGENSOLVE(Hfun, Km + kpt);
    Hm = diag(0);  

    % Block 2~4: 三个卫星点 Hi = Hfun(Km + kpt + Gm(i,:))
    psi_sat = cell(3,1);
    for i = 1:3
        [ei, psi_sat{i}] = vasplib.EIGENSOLVE(Hfun, kpt + Gm(i,:));
        Hm = blkdiag(Hm, diag(ei));
    end
    
    
    % 现在 Hm 是 4N×4N，前 N 行/列是中心，后续每 N×N 是一个卫星
    % 填充中心 ↔ 卫星 的耦合
    %     for i = 1:3
    %         col_idx = (N*i+1) : (N*(i+1));
    %         % 只保留对角方式耦合： <ψ_sat|Vm|ψ0>
    %         Hm(1:N,    col_idx) = psi_sat{i+1}' * Vm * psi_sat{1};
    %         Hm(1:/N,    col_idx) =  Vm ;
    %         Hm(col_idx, 1:N   ) = Hm(1:N, col_idx)';        % Hermitian 共轭
    %     end
    
    M = 3;  % 假设有4个block
    for i = 1:M
        row_idx = (N*(i-1)+1):(N*i);
        for j = i+1:M
            col_idx = (N*(j-1)+1):(N*j);
            
            % 耦合项 <ψ_sat{i}| Vm |ψ_sat{j}>
            Hij = psi_sat{i}' * Vm * psi_sat{j};
%             Hij =  Vm ;           
            % 写入 Hm(i,j) 块
            Hm(row_idx, col_idx) = Hij;
            Hm(col_idx, row_idx) = Hij';  % Hermitian conjugate
        end
    end
    
    
    % 强制全矩阵 Hermiticity（消除数值误差）
    Hm = (Hm + Hm')/2;
end