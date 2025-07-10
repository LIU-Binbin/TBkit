function BCPCAR = BCP(Ham, tensor_index, klist, optionsParallel,options)
arguments
    Ham TBkit
    tensor_index (1,2) double = [2,3];
    klist double = [0,0,0];
    optionsParallel.ncore = 4
    options.selectbands = [];
end
% optionscell = namedargs2cell(options);
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
Nbands = Ham.Basis_num();
if isempty(options.selectbands)
    selectbands = 1:(Nbands/2);
else
    selectbands = options.selectbands;
end
nkpts = size(klist, 1);
BCPCAR = zeros(nkpts,length(selectbands));


% volume = dot(cross(Ham.Rm(1,:),Ham.Rm(2,:)),Ham.Rm(3,:));


%% 批处理计算 BCP_ for 1e9 1000*1000*1000
% batch_size = min(options.batch_size, nkpts);
% nbatch = ceil(nkpts / batch_size);
% fprintf('Total k-points: %d, batch size: %d, total batches: %d\n', nkpts, batch_size, nbatch);
% 
% pb = CmdLineProgressBar('Calculating NCTE: ');  % Progress bar for visualization
% 
tic
% for ibatch = 1:nbatch
%     idx_start = (ibatch - 1)*batch_size + 1;
%     idx_end = min(ibatch*batch_size, nkpts);
%     batch_klist = klist(idx_start:idx_end, :);
%     kpts_this_batch = size(batch_klist, 1);
%     %alpha_mu_batch = zeros(kpts_this_batch, nmu);
%     alpha_mu_batch = zeros(1, nmu);
%     pb.print(ibatch,nbatch);
    if use_parallel
        parfor ki = 1:nkpts
            BCPCAR(ki,:)  =  BCP_k(Ham, tensor_index, klist(ki,:), selectbands);
        end
    else
        for ki = 1:nkpts
            BCPCAR(ki,:)  =   BCP_k(Ham, tensor_index, klist(ki,:), selectbands);
        end
    end

toc


% pb.delete();  % Delete progress bar when done
%% 自动关闭 pool（可选）
if use_parallel
    try delete(gcp('nocreate')); catch, end
end

end