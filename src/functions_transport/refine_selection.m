function [selected_idx, weight, coverage_real] = refine_selection(A_list, optionsAdapt)
% REFINE_SELECTION
% ------------------------------------------------------------------------
% 基于 coarse 阶段计算出的 A(k) 做选点（adapted-k 的核心步骤）
%
% 输入：
%   A_list{k}  : coarse 阶段 A(k) 结果，任意矩阵
%
%   optionsAdapt.method ∈ {"maxabs","frobenius","norm2","hybrid","custom"}
%       maxabs    : max(abs(A(:)))        — 非常适合尖峰（Weyl）
%       frobenius : sum(abs(A(:)))        — 平滑 integrand 更稳定
%       norm2     : 2-norm (largest sv)   — 平滑但可抓局域大特征
%       hybrid    : 0.7*maxabs + 0.3*frob — 推荐最稳定
%       custom    : 使用自定义权重函数 custom_weight_fun(A)
%
%   optionsAdapt.adapt_scaling_factor  : Fortran 风格缩放系数
%   optionsAdapt.coverage_target       : 覆盖率阈值（0.9999 → 99.99%）
%   optionsAdapt.min_keep              : 至少保留多少个点
%
% 输出：
%   selected_idx : 被选中的 k 点下标（列向量）
%   weight       : 每点权重（经过缩放+log 后的最终权重）
%   coverage_real: 选中点实际覆盖率（百分比）
% ------------------------------------------------------------------------

arguments
    A_list cell                       % coarse stage output
    optionsAdapt.method (1,1) string  {mustBeMember(optionsAdapt.method, ...
        ["maxabs","frobenius","norm2","hybrid","custom"])} = "maxabs"
    optionsAdapt.AdapteEnable (1,1) logical = true
    optionsAdapt.adapt_scaling_factor (1,1) double {mustBePositive} = 1.0
    optionsAdapt.coverage_target      (1,1) double {mustBePositive} = 0.9999
    optionsAdapt.min_keep             (1,1) double {mustBePositive} = 32

    % for custom method:
    optionsAdapt.custom_weight_fun = []   % @(A) scalar
end

Nk = numel(A_list);
weight_raw = zeros(Nk,1);

%% ------------------------------------------------------------------------
% (1) 根据策略计算 raw weight
%% ------------------------------------------------------------------------
switch optionsAdapt.method
    case "maxabs"
        for k = 1:Nk
            Ak = A_list{k};
            weight_raw(k) = max(abs(Ak(:)));
        end

    case "frobenius"
        for k = 1:Nk
            Ak = A_list{k};
            weight_raw(k) = sum(abs(Ak(:)));
        end

    case "norm2"
        for k = 1:Nk
            Ak = A_list{k};
            weight_raw(k) = norm(Ak,2);   % largest singular value
        end

    case "hybrid"
        for k = 1:Nk
            Ak = A_list{k};
            w1 = max(abs(Ak(:)));
            w2 = sum(abs(Ak(:)));
            weight_raw(k) = 0.7*w1 + 0.3*w2;
        end

    case "custom"
        if isempty(optionsAdapt.custom_weight_fun)
            error("custom_weight_fun must be provided when method='custom'.");
        end
        f = optionsAdapt.custom_weight_fun;
        for k = 1:Nk
            weight_raw(k) = f(A_list{k});
        end
end

%% ------------------------------------------------------------------------
% 特殊情况：全部为 0
%% ------------------------------------------------------------------------
if max(weight_raw)==0
    selected_idx = 1:Nk;% 全选
    weight = weight_raw;
    coverage_real = 0;
    return;
end

%% ------------------------------------------------------------------------
% (2) Fortran 习惯：scale + log(1 + w)
%% ------------------------------------------------------------------------
asf = optionsAdapt.adapt_scaling_factor;
weight = weight_raw * (asf / max(weight_raw));
weight = log(1 + weight);

total_w = sum(weight);

%% ------------------------------------------------------------------------
% (3) 覆盖率选点（核心）
%% ------------------------------------------------------------------------
[sorted_w, ord] = sort(weight, 'descend');
cum_w = cumsum(sorted_w);

idx_cut = find(cum_w >= optionsAdapt.coverage_target * total_w, 1, 'first');
if isempty(idx_cut)
    idx_cut = Nk;
end

idx_cut = max(idx_cut, round(optionsAdapt.min_keep));

selected = false(Nk,1);
selected(ord(1:idx_cut)) = true;

selected_idx = find(selected);

coverage_real = sum(weight(selected)) / total_w * 100;

end
