function chi_abc_mu = SOAHC_int(Ham, tensor_index, klist, mu_list, options)
% intrinsic 2nd order Anomalous Hall effect
% ref: 10.1103/PhysRevLett.127.277202
% in SI unit, Ampere*Volt^-2*meter for 2D case
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
        Ham = Ham.get_dH_dk_func;
    case "HR"
        Ham = Ham.tjmti_gen();
end
%%
chi_abc_mu = 0;

nkpts = size(klist, 1);
tic
if options.ncore == 1
    for ki = 1:nkpts
        chi_abc_mu = chi_abc_mu + SOAHC_int_k(Ham, tensor_index, klist(ki,:), mu_list, 'eps', eps, 'T', T);
    end
else
    pool = parpool(options.ncore);    
    parfor ki = 1:nkpts
        chi_abc_mu = chi_abc_mu + SOAHC_int_k(Ham, tensor_index, klist(ki,:), mu_list, 'eps', eps, 'T', T);
    end    
    delete(pool)
end
toc
%%
volume = dot(cross(Ham.Rm(1,:),Ham.Rm(2,:)),Ham.Rm(3,:));
chi_abc_mu = chi_abc_mu * constants.charge_C / constants.hbar_eV_s / nkpts / volume;
end