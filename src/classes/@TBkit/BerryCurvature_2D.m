function [BCCAR, Grid] = BerryCurvature_2D(Ham, tensor_index, optionsK, options)
arguments
    Ham TBkit;
    tensor_index (1,2) double = [1 2] % Omega_xy

    optionsK.kstart(1,3) double = [0 0 0]
    optionsK.kdir1 (1,3) double = [1 0 0]
    optionsK.kdir2 (1,3) double = [0 1 0]
    optionsK.Nk1 double = 51
    optionsK.Nk2 double = 51

    options.BAND_index = []; 
    % options.ncore = 1
end
%%
optionsKcell = namedargs2cell(optionsK);
[klist, ~, ~, ~, ~, Grid] = kmeshgen(Ham.Rm, optionsKcell{:}, 'dimension', 2);
%%
switch class(Ham)
    case {"HK", "Htrig"}
    case "HR"
        Ham = Ham.tjmti_gen();
        Ham.Basis_num = Ham.WAN_NUM;
end
% %% start matlab parallel pool (if needed)
% use_parallel = (options.ncore > 1);
% if use_parallel
%     pool = gcp('nocreate');
%     if isempty(pool) || pool.NumWorkers ~= optionsParallel.ncore
%         if ~isempty(pool)
%             delete(pool);
%         end
%         pool = parpool(optionsParallel.ncore);
%     end
% end
%%
mu = 0;
nkpts = size(klist, 1);
BCCAR = zeros(nkpts,1);
%%
clear fft
tic
if isempty(options.BAND_index)
    for ki = 1:nkpts
        BCCAR(ki) = BerryCurvature_k(Ham, tensor_index, klist(ki,:), mu);
    end
else
    for ki = 1:nkpts
        Omega = BerryCurvature_nk(Ham, tensor_index, klist(ki,:));
        BCCAR(ki) = sum(Omega(options.BAND_index));
    end
end
toc
%%
BCCAR = reshape(BCCAR, [optionsK.Nk1, optionsK.Nk2]);
end