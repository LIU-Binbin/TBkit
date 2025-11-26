function refined_klist = refined_klist_gen(klist, selected_idx, mesh_shape, refine_factor, Rm)
% refined_klist_gen
% ---------------------------------------------------------
% 输入：
%   klist         ：原始 coarse 网格
%   selected_idx  ：重要 k 点下标
%   mesh_shape    ：[Nk1 Nk2 Nk3]
%   refine_factor ：整数 refinement 倍数
%   Rm            ：格子矩阵（可选）
%
% 输出：
%   refined_klist ：扩张后的 klist
% ---------------------------------------------------------

if refine_factor == 1
    refined_klist = klist(selected_idx,:);
    return;
end

if isempty(mesh_shape)
    error("refined_klist_gen:MeshShapeMissing", ...
        "mesh_shape must be provided for refinement.");
end
if prod(mesh_shape) ~= size(klist,1)
    error("refined_klist_gen:MeshShapeMismatch", ...
        "product(mesh_shape) != Nk");
end

% ---------------------------------------------------------
% 生成 offsets
% ---------------------------------------------------------
if ~isempty(Rm)
    offsets = refine_recip(mesh_shape, refine_factor, Rm);
else
    step = infer_step(klist, mesh_shape);
    offsets = refine_steps(mesh_shape, refine_factor, step);
end

nref = size(offsets,1);
Ns = numel(selected_idx);

refined_klist = zeros(Ns*nref, size(klist,2));

ptr = 1;
for i = 1:Ns
    base = klist(selected_idx(i),:);
    refined_klist(ptr:ptr+nref-1,:) = base + offsets;
    ptr = ptr + nref;
end
end


%% =============================================================
% 子函数：直接使用你原来的代码
%% =============================================================
function offsets = refine_recip(mesh_shape, factor, Rm)
Nk1=mesh_shape(1); Nk2=mesh_shape(2); Nk3=mesh_shape(3);

sub = ones(1,3);
sub(mesh_shape>1) = factor;

[g1,g2,g3] = ndgrid(0:sub(1)-1,0:sub(2)-1,0:sub(3)-1);

frac = [g1(:)/(sub(1)*Nk1), ...
        g2(:)/(sub(2)*Nk2), ...
        g3(:)/(sub(3)*Nk3)];

Gk = 2*pi * eye(3) / Rm;
offsets = frac * Gk;
end

function offsets = refine_steps(mesh_shape, factor, step)
sub = ones(1,3);
sub(mesh_shape>1) = factor;
[g1,g2,g3] = ndgrid(0:sub(1)-1,0:sub(2)-1,0:sub(3)-1);
offsets = ...
    g1(:).*(step(1,:)/sub(1)) + ...
    g2(:).*(step(2,:)/sub(2)) + ...
    g3(:).*(step(3,:)/sub(3));
end

function step = infer_step(klist, mesh_shape)
Nk1 = mesh_shape(1); Nk2 = mesh_shape(2); Nk3 = mesh_shape(3);
step = zeros(3,size(klist,2));
if Nk1>1, step(1,:) = klist(2,:) - klist(1,:); end
if Nk2>1, step(2,:) = klist(Nk1+1,:) - klist(1); end
if Nk3>1, step(3,:) = klist(Nk1*Nk2+1,:) - klist(1); end
end
