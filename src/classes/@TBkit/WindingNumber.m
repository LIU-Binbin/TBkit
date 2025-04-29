function [WindingNumber,WL] = WindingNumber(Ham_obj,GammaOper,kloop,options)
% WINDINGNUMBER Calculate winding number for topological systems
%
% Syntax:
%   [WindingNumber,WL] = WindingNumber(Ham_obj,GammaOper,kloop)
%   [WindingNumber,WL] = WindingNumber(...,options)
%
% Description:
%   Computes the winding number topological invariant for systems with chiral
%   symmetry, following the formula:
%   N = (1/4πi)∮Tr[Γ·H⁻¹·∇H]·dq
%   where Γ is the chiral symmetry operator.
%
% Input Arguments:
%   Ham_obj - Hamiltonian object (HR, Htrig, HK, or symbolic)
%   GammaOper - Chiral symmetry operator matrix
%   kloop - k-point loop for integration
%   options - Structure containing:
%       plot - Plotting flag (default: false)
%       dir - Directions for calculation (default: [1,2,3])
%
% Output Arguments:
%   WindingNumber - Calculated winding number (scalar)
%   WL - Winding number contributions along loop (array)
%
% Example:
%   Gamma = [0,1;1,0]; kloop = linspace(0,1,100)';
%   [wn,wl] = WindingNumber(Hk,Gamma,kloop);
%
% Reference:
%   https://journals.aps.org/prb/abstract/10.1103/PhysRevB.99.121106
%
% See also: WilsonLoop, TopologicalInvariant
arguments
    Ham_obj;
    GammaOper =[]; % the chiral symmetry oper of Ham_obj
    kloop double = [];
    options.plot logical = false;
    options.dir = [1,2,3];
end
dir = options.dir ;
% prepare dH & dH_dk
switch class(Ham_obj)
    case "HR"
        switch Ham_obj.Type
            case {'sparse'}
                HnumList = reshape(full(cell2mat(Ham_obj.HnumL)),Ham_obj.WAN_NUM,Ham_obj.WAN_NUM,Ham_obj.NRPTS);
            case {'mat'}
                HnumList = Ham_obj.HnumL;
            case {'list'}
                Ham_obj = Ham_obj.rewind();
                HnumList = Ham_obj.HnumL;
        end
        NRPTS_ = Ham_obj.NRPTS;
        vectorList = double(Ham_obj.vectorL);
        vectorList_r = vectorList * Ham_obj.Rm;
        % partial R
        HnumLpx = 1i*pagemtimes(reshape(vectorList_r(:,dir(1)),[1 1 NRPTS_]),HnumList);
        HnumLpy = 1i*pagemtimes(reshape(vectorList_r(:,dir(2)),[1 1 NRPTS_]),HnumList);
        HnumLpz = 1i*pagemtimes(reshape(vectorList_r(:,dir(3)),[1 1 NRPTS_]),HnumList);
        % partial titj  % we dont consider tij mat here check!
        Ham_obj = Ham_obj.tjmti_gen();
        tji_mat_frac = Ham_obj.tjmti{2};
        HnumLpA_tji = tji_mat_frac(:,:,dir(1));
        HnumLpB_tji = tji_mat_frac(:,:,dir(2));
        HnumLpC_tji = tji_mat_frac(:,:,dir(3));
        kloop_cart = kloop * Ham_obj.Gk;
    case {"Htrig"}
        [dH_dkx_fun,dH_dky_fun,dH_dkz_fun] = Ham_diff(Ham_obj);
        HfunTmp = Ham_obj.Hfun;
        kloop_cart = kloop;
    case 'HK'
        syms k_x  k_y k_z real;
        TargetH =  Ham_obj.Hk_num ;
        TargetH_Pk_x = diff(TargetH,k_x);
        TargetH_Pk_y = diff(TargetH,k_y);
        TargetH_Pk_z = diff(TargetH,k_z);
        HfunTmp = matlabFunction(TargetH,'Vars',[k_x k_y k_z]);
        dH_dkx_fun = matlabFunction(TargetH_Pk_x,'Vars',[k_x k_y k_z]);
        dH_dky_fun = matlabFunction(TargetH_Pk_y,'Vars',[k_x k_y k_z]);
        dH_dkz_fun = matlabFunction(TargetH_Pk_z,'Vars',[k_x k_y k_z]);
        kloop_cart = kloop;
    otherwise
        syms k_x k_y k_z real;
        TargetH =  Ham_obj.Hsym ;
        TargetH_Pk_x = diff(TargetH,k_x);
        TargetH_Pk_y = diff(TargetH,k_y);
        TargetH_Pk_z = diff(TargetH,k_z);
        HfunTmp = matlabFunction(TargetH,'Vars',[k_x k_y k_z]);
        dH_dkx_fun = matlabFunction(TargetH_Pk_x,'Vars',[k_x k_y k_z]);
        dH_dky_fun = matlabFunction(TargetH_Pk_y,'Vars',[k_x k_y k_z]);
        dH_dkz_fun = matlabFunction(TargetH_Pk_z,'Vars',[k_x k_y k_z]);
        kloop_cart = kloop;
end
% diff k
% check 1 == end
if ~Oper.isclose(kloop_cart(1,:),kloop_cart(end,:))
    error('loop enforced!');
end
dkloop = (diff(kloop_cart,1,1) + diff(kloop_cart([end-1,1:end-1],:),1,1))/2;
AL = zeros(size(dkloop));
switch class(Ham_obj)
    case "HR"
        for kn = 1:size(dkloop,1)
            % efactor R
            FactorListki = exp(1i*2*pi*vectorList*kloop(kn,:).');
            % HRmat
            HRmat = sum(pagemtimes(HnumList,reshape(FactorListki,[1 1 NRPTS_])),3);
            % pHRmat
            HRmatpA = sum(pagemtimes(HnumLpx,reshape(FactorListki,[1 1 NRPTS_])),3);
            HRmatpB = sum(pagemtimes(HnumLpy,reshape(FactorListki,[1 1 NRPTS_])),3);
            HRmatpC = sum(pagemtimes(HnumLpz,reshape(FactorListki,[1 1 NRPTS_])),3);
            % efactor orb
            kjiL_A =  tji_mat_frac(:,:,1).*kloop(kn,1);
            kjiL_B =  tji_mat_frac(:,:,2).*kloop(kn,2);
            kjiL_C =  tji_mat_frac(:,:,3).*kloop(kn,3);
            Hmat_tji = exp(1i*2*pi*(kjiL_A+kjiL_B+kjiL_C));
            %
            Hmat_tjipA = Hmat_tji.* HnumLpA_tji;
            Hmat_tjipB = Hmat_tji.* HnumLpB_tji;
            Hmat_tjipC = Hmat_tji.* HnumLpC_tji;
            %
            vxk = HRmatpA.*Hmat_tji + HRmat.*Hmat_tjipA;% vx partial_A_tmp
            vyk = HRmatpB.*Hmat_tji + HRmat.*Hmat_tjipB;% vy partial_B_tmp
            vzk = HRmatpC.*Hmat_tji + HRmat.*Hmat_tjipC;% vx partial_A_tmp
            H = HRmat.*Hmat_tji;
            Ax = trace(GammaOper/H*vxk);
            Ay = trace(GammaOper/H*vyk);
            Az = trace(GammaOper/H*vzk);
            AL(kn,:)= [Ax,Ay,Az];
        end
    otherwise
        for kn = 1:size(dkloop,1)
            kx = kloop_cart(kn,1); ky = kloop_cart(kn,2); kz = kloop_cart(kn,3);
            dH_dkx = dH_dkx_fun(kx,ky,kz);dH_dky = dH_dky_fun(kx,ky,kz);
            dH_dkz = dH_dkz_fun(kx,ky,kz);H = HfunTmp(kx,ky,kz);
            Ax = trace(GammaOper/H*dH_dkx);
            Ay = trace(GammaOper/H*dH_dky);
            Az = trace(GammaOper/H*dH_dkz);
            AL(kn,:)= [Ax,Ay,Az];
        end
end
WL = dot(AL,dkloop,2)/(4i*pi);
WindingNumber = sum(WL);
end