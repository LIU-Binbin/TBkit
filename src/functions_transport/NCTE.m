function alpha_mu = NCTE(Ham, tensor_index, klist, mu_list, optionsParallel, options)
% PARAMETERS
% ----------
% Ham    : TBkit HK/HR object
% tensor_index : [1×3] or [] for all 27 components
% klist  : [Nk × 3] k-point list
% mu_list: [1 × Nμ] chemical potentials
%
% optionsParallel.ncore : number of cores
%
% options.T    : temperature (scalar)
% options.eps  : degeneracy threshold
% options.batch_size : (reserved for future multi-batch integration)
%
% RETURNS
% -------
% alpha_mu :
%     If tensor_index is single (a,b,c): [1 × Nμ]
%     If tensor_index = []:             [27 × Nμ]
%
% ---------------------------------------------------------------
arguments
    Ham TBkit
    tensor_index  double = [];
    klist double = [];
    mu_list double = [];
    optionsParallel.ncore = 4
    options.T = 50 % Kelvin
    options.eps = 1e-4
    options.batch_size = 1e6  % 默认批次大小 100*100*100
    options.Coeffs = [1 1 -1];
    options.adapt.enable (1,1) logical = false
    options.adapt.importance_factor (1,1) double {mustBePositive} = 0.2
    options.adapt.min_keep (1,1) double {mustBePositive} = 32
    options.adapt.reuse_samples (1,1) logical = true
    options.adapt.mesh_shape double = []
    options.adapt.refine_factor (1,1) double {mustBePositive} = 1
end

% optionscell = namedargs2cell(options);
%% prepare dH_dk
switch class(Ham)
    case "HK"
    case "HR"
        Ham = Ham.tjmti_gen();
end

clear fft;
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
% ------------- Basic info -----------------------
Nmu = numel(mu_list);
Nk  = size(klist,1);

% Output shape:
if isempty(tensor_index)
    alpha_mu = zeros(27, Nmu);
else
    alpha_mu = zeros(1, Nmu);
end

T   = options.T;
eps = options.eps;

% -------- Volume factor and prefactor ----------
volume = dot(cross(Ham.Rm(1,:),Ham.Rm(2,:)), Ham.Rm(3,:));

const_factor = (constants.charge_C / T / constants.hbar_eV_s) ...
    * (volume / Nk) ...
    / constants.hbar_eV_s;
%%
% ------------------------------------------------
%   MAIN LOOP OVER k-POINTS
% ------------------------------------------------
k_point_fun = @(kpt) NCTE_k(Ham, tensor_index, kpt, mu_list, T, eps, options.Coeffs);

alpha_mu = kloop_accumulate(klist, k_point_fun, alpha_mu, optionsParallel, options.adapt);

%% -------- Multiply physical prefactor ----------
alpha_mu = alpha_mu * const_factor;
end
