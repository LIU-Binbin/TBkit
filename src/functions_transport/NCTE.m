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
% optionscell = namedargs2cell(options);
%% prepare dH_dk
switch class(Ham)
    case "HK"
        
    case "HR"
        Ham = Ham.tjmti_gen();
end
%% start matlab parallel pool (if needed)
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
const_factor = (constants.charge_C / T / constants.hbar_eV_s) * (volume / nkpts) / constants.hbar_eV_s;
%%
if use_parallel
    parfor ki = 1:nkpts
        alpha_mu  = alpha_mu + NCTE_k(Ham, tensor_index, klist(ki,:), mu_list, T, eps);
    end
    %delete(pool)
else
    for ki = 1:nkpts
        alpha_mu  = alpha_mu + NCTE_k(Ham, tensor_index, klist(ki,:), mu_list, T, eps);
    end
end
% toc
%%
alpha_mu = alpha_mu * const_factor;
end