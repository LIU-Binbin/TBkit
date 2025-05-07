function GREENCAR = GREENCAR_gen(w_list,eta,H00_H01_cell_list_1,H00_H01_cell_list_2,mode,mu_max)
%GREENCAR_GEN Generate Green's function data for multiple k-points
%
%   Syntax:
%       GREENCAR = GREENCAR_gen(w_list,eta,H00_H01_cell_list_1,H00_H01_cell_list_2,mode,mu_max)
%
%   Description:
%       Computes Green's functions using different algorithms:
%       - Recursive Green's function (Green_iter)
%       - T-matrix methods (Tmatrix_iter, Tmatrix_inv)
%
%   Inputs:
%       w_list              - Frequency/energy points
%       eta                 - Broadening parameter
%       H00_H01_cell_list_1 - Hamiltonian diagonal blocks
%       H00_H01_cell_list_2 - Hamiltonian off-diagonal blocks
%       mode                - Calculation method:
%                             'Green_iter', 'Tmatrix_iter', 'Tmatrix_inv'
%       mu_max              - Maximum recursion depth (default=100)
%
%   Output:
%       GREENCAR - Structure containing bulk and surface Green's functions
%
%   See also: GW_iter, Tmatrix_gen
arguments
    w_list 
    eta 
    H00_H01_cell_list_1 
    H00_H01_cell_list_2 
    mode = 'Green_iter';
    mu_max =  100;
end

switch mode
    case 'Green_iter'
        %H00_H01_cell_list = H00_H01_cell_list;
        WAN_NUM = size(H00_H01_cell_list_1,1);
        kn =size( H00_H01_cell_list_2,3);
        Nsize1 = length(w_list);
        if Nsize1 == 1
            % arc
        end
        Nsize2 = kn;
        GREENCAR1{Nsize1,Nsize2} = zeros(WAN_NUM) ;
        GREENCAR2{Nsize1,Nsize2} = zeros(WAN_NUM) ;
        GREENCAR3{Nsize1,Nsize2} = zeros(WAN_NUM) ;
        count = 0;
        pb = CmdLineProgressBar('Green Solving ');
        Nsize = Nsize1*Nsize2;
        for i = 1:Nsize1
            for j= 1:Nsize2
                count = count +1;
                pb.print(count,Nsize);
                w =w_list(i);
                H00 = H00_H01_cell_list_1(:,:,j);
                H01 = H00_H01_cell_list_2(:,:,j);
                [GREENCAR3{i,j},GREENCAR1{i,j},GREENCAR2{i,j}] = HR.GW_iter(H00,H01,w,eta,mu_max);
            end
        end
        pb.delete()
        GREENCAR.bulk = GREENCAR1;
        GREENCAR.surf_l = GREENCAR2;
        GREENCAR.surf_r = GREENCAR3;
    case 'Tmatrix_iter'
        %H00_H01_cell_list = H00_H01_cell_list;
        WAN_NUM = size(H00_H01_cell_list_1,1);
        kn =size( H00_H01_cell_list_2,3);
        Nsize1 = length(w_list);
        Nsize2 = kn;
        GREENCAR1{Nsize1,Nsize2} = zeros(WAN_NUM) ;
        GREENCAR2{Nsize1,Nsize2} = zeros(WAN_NUM) ;
        GREENCAR3{Nsize1,Nsize2} = zeros(WAN_NUM) ;
        for i = 1:Nsize1
            for j= 1:Nsize2
                w =w_list(i);
                H00 = H00_H01_cell_list_1(:,:,j);
                H01 = H00_H01_cell_list_2(:,:,j);
                [GREENCAR3{i,j},GREENCAR1{i,j},GREENCAR2{i,j}] = HR.Tmatrix_iter(H00,H01,w,eta);
                % [GREENCAR1{i,j}, GREENCAR2{i,j},GREENCAR3{i,j}]=   Tmatrix2Green00(Tmatrix,H00,H01,w,eta);
            end
        end
        GREENCAR.bulk = GREENCAR1;
        GREENCAR.surf_l = GREENCAR2;
        GREENCAR.surf_r = GREENCAR3;
        %                 case  't'
        %                     if iscell(H00_H01_cell_list)
        %                         disp('do nothing');
        %                     elseif isstruct(H00_H01_cell_list)
        %                         TBkitobj = H00_H01_cell_list;
        %                         H00_H01_cell_list =H_cell_list_gen(TBkitobj);
        %                     end
        %                     Nsize1 = length(w_list);
        %                     Nsize2 = length(H00_H01_cell_list);
        %                     %w_list= linspace(w_range(1),w_range(2),w_number);
        %                     WAN_NUM = length(H00_H01_cell_list{1});
        %                     GREENCAR{Nsize1,Nsize2} = zeros(WAN_NUM) ;
        %                     for i = 1:Nsize1
        %                         for j= 1:Nsize2
        %                             GREENCAR{i,j} =  GREENCAR_one_gen(H00_H01_cell_list{j},w_list(i),eta,'TB');
        %                         end
        %                     end
        %                 case 'Tmatrix_eig'
        %                     %H00_H01_cell_list = H00_H01_cell_list;
        %                     WAN_NUM = length(H00_H01_cell_list{1,1});
        %                     [~,kn] =size( H00_H01_cell_list);
        %                     Nsize1 = length(w_list);
        %                     Nsize2 = kn;
        %                     GREENCAR{Nsize1,Nsize2} = zeros(WAN_NUM) ;
        %                     for i = 1:Nsize1
        %                         for j= 1:Nsize2
        %                             w =w_list(i);
        %                             H00 = H00_H01_cell_list{1,j};
        %                             H01 = H00_H01_cell_list{2,j};
        %                             Tmatrix = HR.Tmatrix_gen(H00,H01,w,0);
        %                             GREENCAR{i,j} = Tmatrix2pband(Tmatrix,H00,H01,w);
        %                         end
        %                     end
    case 'Tmatrix_inv'
        %H00_H01_cell_list = H00_H01_cell_list;
        WAN_NUM = size(H00_H01_cell_list_1,1);
        kn =size( H00_H01_cell_list_2,3);
        Nsize1 = length(w_list);
        Nsize2 = kn;
        GREENCAR1{Nsize1,Nsize2} = zeros(WAN_NUM) ;
        GREENCAR2{Nsize1,Nsize2} = zeros(WAN_NUM) ;
        GREENCAR3{Nsize1,Nsize2} = zeros(WAN_NUM) ;
        for i = 1:Nsize1
            for j= 1:Nsize2
                w =w_list(i);
                H00 = H00_H01_cell_list_1(:,:,j);
                H01 = H00_H01_cell_list_2(:,:,j);
                Tmatrix = HR.Tmatrix_gen(H00,H01,w,eta);
                [GREENCAR1{i,j}, GREENCAR2{i,j},GREENCAR3{i,j}]=   HR.Tmatrix2Green00(Tmatrix,H00,H01,w,eta);
            end
        end
        GREENCAR.bulk = GREENCAR1;
        GREENCAR.surf_l = GREENCAR2;
        GREENCAR.surf_r = GREENCAR3;


end
end