function WAVECAR_loop = Topo1DpreWAVECAR(TBkitobj,klist,options)
% TOPO1DPREWAVECAR Prepare WAVECAR for topological calculations in 1D
%
% Syntax:
%   WAVECAR_loop = Topo1DpreWAVECAR(TBkitobj,klist)
%   WAVECAR_loop = Topo1DpreWAVECAR(TBkitobj,klist,options)
%
% Description:
%   Prepares wavefunction data (WAVECAR) for topological calculations
%   along a 1D k-path for different TBkit object types.
%
% Input Arguments:
%   TBkitobj - Tight-binding object (Htrig, HK, or HR type)
%   klist - List of k-points in reciprocal space
%   options - Optional parameters:
%       BAND_index - Band indices to include (default: half of total bands)
%       ax - Axis handle for plotting
%       plot - Flag to enable plotting (default: false)
%       LWAVE - Flag for wavefunction processing (default: false)
%
% Output Arguments:
%   WAVECAR_loop - Wavefunction data array (Norb x Nband x Nk)
%
% Example:
%   k = linspace(0,1,10)';
%   W = Topo1DpreWAVECAR(Hk,k);
arguments
    TBkitobj
    klist
    options.BAND_index = [];
    options.ax = handle([]);
    options.plot = false;
    options.LWAVE = false;
end
% --- nargin
switch class(TBkitobj)
    case {'Htrig','HK'}
        %
    case 'HR'
        TBkitobj.Basis_num = TBkitobj.WAN_NUM;
end
if isempty(options.BAND_index)
    set_divide = 2;
    switch class(TBkitobj)
        case {'Htrig','HK','HR'}
            BAND_index = 1:(TBkitobj.Basis_num/set_divide);
        case ''
            BAND_index = 1:(TBkitobj.Basis_num/set_divide);
        otherwise
            BAND_index = 1:(TBkitobj.Basis_num/set_divide);
    end
else
    BAND_index = options.BAND_index;
end
switch class(TBkitobj)
    case {'Htrig','HK'}
        klist_tmp = klist;
        [~,WAVECAR_loop] = TBkitobj.EIGENCAR_gen(...
            'klist',klist_tmp,'printmode',false);
        WAVECAR_loop = WAVECAR_loop(:,BAND_index,:);
        % The last bloch state is the same as the first up to a phase factor
        WAVECAR_loop(:,:,end) = WAVECAR_loop(:,:,1).* exp(-1i*(TBkitobj.orbL*TBkitobj.Rm*(klist_tmp(end,:)-klist_tmp(1,:)).'));
    case 'HR'
        klist_tmp = klist;
        [~,WAVECAR_loop] = TBkitobj.EIGENCAR_gen(...
            'klist',klist_tmp,...
            'convention','II','printmode',false);
        % If we use convention II, each wavefactor should
        % give back the factor
        % C^{nk}_j = C^{nk}_j_tilde * e^{-ik·tj}.
        WAVECAR_loop = WAVECAR_loop(:,BAND_index,:);
        % normalize phases to get u instead of phi
        for j =1:size(WAVECAR_loop,3)
            WAVECAR_loop(:,:,j) = WAVECAR_loop(:,:,j).* exp(-2i*pi*(TBkitobj.orbL*klist_tmp(j,:).'));
        end
        % The last bloch state is the same as the first up to a phase factor
        WAVECAR_loop(:,:,end) = WAVECAR_loop(:,:,1).* exp(-2i*pi*(TBkitobj.orbL*(klist_tmp(end,:)-klist_tmp(1,:)).'));
    otherwise
        klist_tmp = klist;
        [~,WAVECAR_loop] = TBkit.EIGENSOLVE(TBkitobj.Hfun,klist_tmp,TBkitobj.Basis_num);
        WAVECAR_loop = WAVECAR_loop(:,BAND_index,:);
        % The last bloch state is the same as the first up to a phase factor
        WAVECAR_loop(:,:,end) = WAVECAR_loop(:,:,1).* exp(-1i*(TBkitobj.orbL*TBkitobj.Rm*(klist_tmp(end,:)-klist_tmp(1,:)).'));
end

end