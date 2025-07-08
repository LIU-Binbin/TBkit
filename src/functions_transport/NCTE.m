function alpha_mu = NCTE(Ham, tensor_index, klist, mu_list, optionsParallel,options)
arguments
    Ham TBkit
    tensor_index (1,3) double
    klist double
    mu_list double
    optionsParallel.ncore = 4
    options.T = 50 % Kelvin
    options.eps = 1e-4
    options.batch_size = 1e6  % 默认批次大小 100*100*100
end
optionscell = namedargs2cell(options);
%% prepare dH_dk
switch class(Ham)
    case "HK"
        
    case "HR"
        Ham = Ham.tjmti_gen();
end
%% 启动并行池（如需要）
use_parallel = (optionsParallel.ncore > 1);
if use_parallel
    pool = gcp('nocreate');
    if isempty(pool) || pool.NumWorkers ~= optionsParallel.ncore
        if ~isempty(pool)
            delete(pool);
        end
        pool = parpool(optionsParallel.ncore);
    end
end
%%
nmu = length(mu_list);
alpha_mu = zeros(1, nmu); 
nkpts = size(klist, 1);
volume = dot(cross(Ham.Rm(1,:),Ham.Rm(2,:)),Ham.Rm(3,:));
T = options.T;
eps = options.eps;
const_factor = constants.charge_C / constants.hbar_eV_s^2 / nkpts / volume / T;

%% 批处理计算 NCTE_k for 1e9 1000*1000*1000
batch_size = min(options.batch_size, nkpts);
nbatch = ceil(nkpts / batch_size);
fprintf('Total k-points: %d, batch size: %d, total batches: %d\n', nkpts, batch_size, nbatch);

pb = CmdLineProgressBar('Calculating NCTE: ');  % Progress bar for visualization

tic
for ibatch = 1:nbatch
    idx_start = (ibatch - 1)*batch_size + 1;
    idx_end = min(ibatch*batch_size, nkpts);
    batch_klist = klist(idx_start:idx_end, :);
    kpts_this_batch = size(batch_klist, 1);
    alpha_mu_batch = zeros(kpts_this_batch, nmu);

    if use_parallel
        parfor ki = 1:kpts_this_batch
            alpha_mu_batch(ki,:) = NCTE_k(Ham, tensor_index, batch_klist(ki,:), mu_list, T, eps);
        end
    else
        for ki = 1:kpts_this_batch
            alpha_mu_batch(ki,:) = NCTE_k(Ham, tensor_index, batch_klist(ki,:), mu_list, T, eps);
        end
    end

    alpha_mu = alpha_mu + sum(alpha_mu_batch, 1);

    % 可选保存中间结果，防止崩溃丢失
    % save(sprintf('alpha_mu_batch_%d.mat', ibatch), 'alpha_mu_batch', 'idx_start', 'idx_end');
end
toc

% tic
% 
%     if use_parallel
%         parfor ki = 1:nkpts
%             alpha_mu = alpha_mu + NCTE_k(Ham, tensor_index, klist(ki,:), mu_list, T, eps);
%             %alpha_mu = alpha_mu + alpha_mu_temp;
%         end
%     else
%         for ki = 1:nkpts
%             alpha_mu = alpha_mu + NCTE_k(Ham, tensor_index, klist(ki,:), mu_list, T, eps);
%             %alpha_mu = alpha_mu + alpha_mu_temp;
%         end
%     end
% 
% toc

% pb.delete();  % Delete progress bar when done
%%
alpha_mu = alpha_mu * const_factor;
%% 自动关闭 pool（可选）
if use_parallel
    try delete(gcp('nocreate')); catch, end
end

end