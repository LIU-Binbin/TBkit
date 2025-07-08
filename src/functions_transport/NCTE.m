function alpha_mu = NCTE(Ham, tensor_index, klist, mu_list, options)
arguments
    Ham TBkit
    tensor_index (1,3) double
    klist double
    mu_list double
    options.ncore = 4
    options.T = 50 % Kelvin
    options.eps = 1e-4
end
eps = options.eps;
T = options.T;
%% prepare dH_dk
switch class(Ham)
    case "HK"
        
    case "HR"
        Ham = Ham.tjmti_gen();
end
%%
alpha_mu = 0;

nkpts = size(klist, 1);

pool = parpool(options.ncore);
tic
% pb = CmdLineProgressBar('Calculating NCTE: ');  % Progress bar for visualization
parfor ki = 1:nkpts
    alpha_mu = alpha_mu + NCTE_k(Ham, tensor_index, klist(ki,:), mu_list, 'eps', eps, 'T', T);
    % pb.print(ki, nkpts, " kpoints done");
end
% pb.delete();  % Delete progress bar when done
toc
delete(pool)
%%
volume = dot(cross(Ham.Rm(1,:),Ham.Rm(2,:)),Ham.Rm(3,:));
alpha_mu = alpha_mu * constants.charge_C / constants.hbar_eV_s^2 / nkpts / volume / T;
end