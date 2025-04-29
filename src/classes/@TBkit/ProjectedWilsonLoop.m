function [BFCAR,WEIGHTCAR,klist_l] = ProjectedWilsonLoop(TBkitobj,options)
%PROJECTEDWILSONLOOP Compute Wilson loop with orbital projections
%
%   Syntax:
%       [BFCAR,WEIGHTCAR,klist_l] = ProjectedWilsonLoop(TBkitobj,options)
%
%   Description:
%       Calculates Wilson loop spectrum with optional symmetry projections.
%       Supports both Cartesian and reciprocal coordinates.
%
%   Inputs:
%       TBkitobj - Tight-binding model object
%       options  - Calculation options structure:
%                  BAND_index, knum_int, knum_evol, etc.
%
%   Outputs:
%       BFCAR      - Wilson loop eigenvalues
%       WEIGHTCAR  - Projection weights
%       klist_l    - k-path coordinates
%
%   See also: WilsonLoop, nested_WilsonLoop
arguments
    TBkitobj ;
    options.BAND_index = [];
    options.knum_int    = 31;
    options.knum_evol   = 51;
    options.kstart      = [0,0,0];
    options.kintegral   = [0,1,0];
    options.kevolution  = [1,0,0];
    options.cartesian = false;
    options.dir_seq = [1,2,3];
    options.dir_start = 'k_z';
    options.fig = false;
    options.plot = false;
    options.Oper = [];
    options.ProjectionMethod = 'sign';
    options.ProjectionStruct = struct('field','imag');
end
% --- nargin
if isempty(options.Oper)
    project = false;
    set_divide = 2;
else
    project = true;
    set_divide = 4;
    OperObj = options.Oper; % only support one Oper now
    OperU = roundn(OperObj.U,-6);
    [AOperU,UOperU] = eig(OperU);
    switch options.ProjectionMethod
        case 'sign'
            if strcmp(options.ProjectionStruct.field,'imag')
                [AOperU,Ucheck] = park.sorteig(imag(UOperU),AOperU);
                %ProjectionFunction = @(WAVECAR) sign(imag(OperObj.EIGEN(WAVECAR)));
            elseif strcmp(options.ProjectionStruct.field,'real')
                [AOperU,Ucheck] = park.sorteig(real(UOperU),AOperU);
                %ProjectionFunction = @(WAVECAR) sign(real(OperObj.EIGEN(WAVECAR)));
            end
        otherwise

    end
end
if isempty(options.BAND_index)
    switch class(TBkitobj)
        case {'Htrig','HK'}

            BAND_index = 1:TBkitobj.Basis_num/set_divide;
        case 'HR'
            BAND_index = 1:TBkitobj.WAN_NUM/set_divide;
    end
else
    BAND_index = options.BAND_index;
end
%
if options.cartesian
    if options.plot
        [fig,ax] = TBkit.BZplot(TBkitobj.Rm,'color','r','alpha',0.1);
        view(ax,1.465409092282576e+02,16.15644804763928);
    end
    Gk_ = TBkit.CartisianMat(TBkitobj.Gk,options.dir_seq,options.dir_start);
    if options.plot
        [fig,ax] = TBkit.BZplot(Gk_,'Gk',true,'color','y','alpha',0.3,'ax',ax,'fig',fig);
        view(ax,1.465409092282576e+02,16.15644804763928);
        disp(Gk_);
    end
    kstart_s  = options.kstart * Gk_ /TBkitobj.Gk;
else
    Gk_ = TBkitobj.Gk;
    kstart_s = options.kstart;
end
%
[klist_r_1,klist_s_1,~,~] =...
    TBkit.kpathgen([[0,0,0];options.kevolution],options.knum_evol,Gk_,TBkitobj.Gk);
klist_l = zeros(size(klist_s_1,1),1);
%             klist_l(1) = sum(sign(klist_s_1(1,:)))*norm(klist_s_1(1,:)*(eye(3)*2*pi));
klist_s_1_ = klist_s_1 ;
normklist_l = norm(options.kevolution)/norm(klist_s_1(end,:));
for i = 1:size(klist_s_1_,1)
    klist_l(i) = norm(klist_s_1_(i,:)*(eye(3)*2*pi))*normklist_l;
end
klist_l = klist_l + sum(sign(kstart_s))* norm(kstart_s*(eye(3)*2*pi))*normklist_l;
[klist_r_2,klist_s_2,~,~] =...
    TBkit.kpathgen([[0,0,0];options.kintegral],options.knum_int,Gk_,TBkitobj.Gk);
BFCAR = zeros(length(BAND_index)*2,options.knum_evol);
WEIGHTCAR = BFCAR;
BF_WAVECAR = zeros(length(BAND_index),length(BAND_index),options.knum_evol);
kstart_r = options.kstart*Gk_;
switch class(TBkitobj)
    case {'Htrig','HK'}
        for i = 1:options.knum_evol
            klist_tmp = klist_r_1(i,:)+klist_r_2+kstart_r;
            [~,WAVECAR_loop1] = TBkitobj.EIGENCAR_gen(...
                'klist',klist_tmp,'show',options.plot);
            WAVECAR_loop_tmp = WAVECAR_loop1(:,BAND_index,:);
            if project
                WEIGHTCAR_SYM = ProjectionFunction(WAVECAR_loop_tmp);
                kn = size(WAVECAR_loop_tmp,3);
                Norb = size(WAVECAR_loop_tmp,1);
                switch options.ProjectionMethod
                    case 'sign'
                        LABELCAR1 = logical(WEIGHTCAR_SYM > 0);
                        LABELCAR2 = logocal(WEIGHTCAR_SYM < 0);
                        PlusNum = sum(LABELCAR1(:,round(end/2)));
                        MinusNum = sum(LABELCAR2(:,round(end/2)));
                        WAVECAR_loop_tmp1 = zeros(Norb,PlusNum,kn);
                        WAVECAR_loop_tmp2 = zeros(Norb,MinusNum,kn);
                        for j = 1:kn
                            WAVECAR_loop_tmp1 = WAVECAR_loop_tmp(:,LABELCAR1(:,j),j);
                            WAVECAR_loop_tmp2 = WAVECAR_loop_tmp(:,LABELCAR2(:,j),j);
                        end
                        [BFCAR1,~] = TBkit.wancenter_1D(WAVECAR_loop_tmp1);% [BFCAR1,BFWAVECAR1]
                        [BFCAR2,~] = TBkit.wancenter_1D(WAVECAR_loop_tmp2);% [BFCAR1,BFWAVECAR2]
                        BFCAR(:,i)=[BFCAR1;BFCAR2];
                        WEIGHTCAR(:,i) = [ones(PlusNum,1);-ones(MinusNum,1)];
                    otherwise
                end
            else
                [BFCAR(:,i),BF_WAVECAR(:,:,i)] = TBkit.wancenter_1D(WAVECAR_loop_tmp);
            end
        end
    case 'HR'
        LABELCAR1 = 1:TBkitobj.WAN_NUM/2;
        LABELCAR2 = (TBkitobj.WAN_NUM/2+1):TBkitobj.WAN_NUM;
        for i = 1:options.knum_evol
            klist_tmp = klist_s_1(i,:)+klist_s_2+kstart_s;
            if project
                if options.plot
                    [~,WAVECAR_loop1,ax] = TBkitobj.EIGENCAR_gen(...
                        'klist',klist_tmp,...
                        'convention','I','show',options.plot,'ax',ax,'Umat',AOperU,'subband',LABELCAR1);
                    drawnow;
                else
                    [~,WAVECAR_loop1] = TBkitobj.EIGENCAR_gen(...
                        'klist',klist_tmp,...
                        'convention','I','printmode',false,'Umat',AOperU,'subband',LABELCAR1);
                end
                %error('debug');
                if options.plot
                    [~,WAVECAR_loop2,ax] = TBkitobj.EIGENCAR_gen(...
                        'klist',klist_tmp,...
                        'convention','I','show',options.plot,'ax',ax,'Umat',AOperU,'subband',LABELCAR2);
                    drawnow;
                else
                    [~,WAVECAR_loop2] = TBkitobj.EIGENCAR_gen(...
                        'klist',klist_tmp,...
                        'convention','I','printmode',false,'Umat',AOperU,'subband',LABELCAR2);
                end
                % If we use convention II, each wavefactor should
                % give back the factor
                % C^{nk}_j = C^{nk}_j_tilde * e^{-ik·tj}.
                % not selected band here?
                WAVECAR_loop_tmp1 = WAVECAR_loop1(:,BAND_index,:);
                WAVECAR_loop_tmp2 = WAVECAR_loop2(:,BAND_index,:);
                % normalize phases to get u instead of phi
                %for j =1:size(WAVECAR_loop_tmp,3)
                %    WAVECAR_loop_tmp(:,:,j) = WAVECAR_loop_tmp(:,:,j).* exp(-2i*pi*(TBkitobj.orbL*klist_tmp(j,:).'));
                %end
                % The last bloch state is the same as the first up to a phase factor
                WAVECAR_loop_tmp1(:,:,end) = WAVECAR_loop_tmp1(:,:,1).* exp(-2i*pi*(TBkitobj.orbL*(klist_tmp(end,:)-klist_tmp(1,:)).'));
                WAVECAR_loop_tmp2(:,:,end) = WAVECAR_loop_tmp2(:,:,1).* exp(-2i*pi*(TBkitobj.orbL*(klist_tmp(end,:)-klist_tmp(1,:)).'));
                %WAVECAR_loop_tmp(:,:,end) = WAVECAR_loop_tmp(:,:,1).* exp(-2i*pi*(TBkitobj.orbL*(klist_tmp(end,:)-klist_tmp(1,:)).'));

                %LABELCAR1 = logical(WEIGHTCAR_SYM > 0);
                %LABELCAR2 = logical(WEIGHTCAR_SYM < 0);
                PlusNum = length(BAND_index);
                MinusNum = length(BAND_index);
                %PlusNum = sum(LABELCAR1(:,round(end/2)));
                %MinusNum = sum(LABELCAR2(:,round(end/2)));
                %for j = 1:kn
                %    WAVECAR_loop_tmp1 = WAVECAR_loop_tmp(:,LABELCAR1(:,j),j);
                %    WAVECAR_loop_tmp2 = WAVECAR_loop_tmp(:,LABELCAR2(:,j),j);
                %end
                WAVECAR_loop_tmp1 = WAVECAR_loop_tmp1(:,BAND_index,:);
                WAVECAR_loop_tmp2 = WAVECAR_loop_tmp2(:,BAND_index,:);
                [BFCAR1,~] = TBkit.wancenter_1D(WAVECAR_loop_tmp1);% [BFCAR1,BFWAVECAR1]
                [BFCAR2,~] = TBkit.wancenter_1D(WAVECAR_loop_tmp2);% [BFCAR1,BFWAVECAR2]
                BFCAR(:,i)=[BFCAR1;BFCAR2];
                WEIGHTCAR(:,i) = [ones(PlusNum,1);-ones(MinusNum,1)];
            else
                klist_tmp = klist_s_1(i,:)+klist_s_2+kstart_s;
                if options.plot
                    [~,WAVECAR_loop,ax] = TBkitobj.EIGENCAR_gen(...
                        'klist',klist_tmp,...
                        'convention','II','show',options.plot,'ax',ax);
                    drawnow;
                else
                    [~,WAVECAR_loop] = TBkitobj.EIGENCAR_gen(...
                        'klist',klist_tmp,...
                        'convention','II','printmode',false);
                end
                % If we use convention II, each wavefactor should
                % give back the factor
                % C^{nk}_j = C^{nk}_j_tilde * e^{-ik·tj}.
                WAVECAR_loop_tmp = WAVECAR_loop(:,BAND_index,:);
                % normalize phases to get u instead of phi
                for j =1:size(WAVECAR_loop_tmp,3)
                    WAVECAR_loop_tmp(:,:,j) = WAVECAR_loop_tmp(:,:,j).* exp(-2i*pi*(TBkitobj.orbL*klist_tmp(j,:).'));
                end
                % The last bloch state is the same as the first up to a phase factor
                WAVECAR_loop_tmp(:,:,end) = WAVECAR_loop_tmp(:,:,1).* exp(-2i*pi*(TBkitobj.orbL*(klist_tmp(end,:)-klist_tmp(1,:)).'));
                [BFCAR(:,i),BF_WAVECAR(:,:,i)] = TBkit.wancenter_1D(WAVECAR_loop_tmp);
            end
        end
end
end