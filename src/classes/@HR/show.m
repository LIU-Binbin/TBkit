function ax = show(H_hr, mode, options)
% Visualizes HR object in different modes: POSCAR, NRPTS, HOPPING
%
% Inputs:
%   H_hr    - HR object containing crystal structure information
%   mode    - Visualization mode {'POSCAR','NRPTS','HOPPING'} (default: 'HOPPING')
%   options - Configuration parameters for visualization
%
% Outputs:
%   fig     - Figure handle
%   ax      - Axes handle

arguments
    H_hr HR; % HR object with structural data
    mode char {mustBeMember(mode,{'POSCAR','NRPTS','HOPPING'})} = 'HOPPING'; % Visualization mode
    options.scale = 1; % Structure scaling factor
    options.atomscale = 0.5; % Atom display scaling
    options.TwoD = false; % 2D projection flag
    options.vectorList = []; % Custom vector list for NRPTS
    options.fast = true; % Fast rendering mode
    options.ax = []; % Existing axes handle
    options.Select = []; % Hopping selection filter
    options.cmap = @parula; % Color mapping function
    options.Title = ''; % Figure title
end
% Prepare axes for plotting
if isempty(options.ax)
    [~,ax]= Figs(1, 1); % Create new figure
else
    if isvalid(options.ax)
        ax = options.ax; % Use provided axis handle
    else
        [~,ax]= Figs(1, 1); % Create new figure 
    end
end

% Scale lattice vectors
Rm_ = H_hr.Rm * options.scale;

% Handle orthogonal vector conversion
if H_hr.vectorhopping
    H_hr = H_hr.GenfromOrth(); % Convert to orthogonal vectors
end

% Dimension validation
if H_hr.Dim > 3
    error('Visualization only supports dimensions <= 3');
end

% Convert orbital positions to double precision
H_hr.orbL = double(H_hr.orbL);

switch mode
    case 'POSCAR' % Basic structure plot
        ax = POSCAR_plot(Rm_, H_hr.sites, H_hr.Atom_name, H_hr.Atom_num,'ax',ax);

    case 'NRPTS' % Repeated unit plot
        % Generate/use vector list
        if isempty(options.vectorList)
            vectorList = double(unique(H_hr.vectorL(:,1:H_hr.Dim), 'rows'));
        else
            vectorList = options.vectorList;
        end

        NRPTS_ = size(vectorList, 1); % Number of repeated units
        positions = H_hr.orbL; % Base positions

        % Build expanded site list
        sites_ = [];
        if isempty(H_hr.elementL)
            elementList = 6 * ones(size(positions, 1), 1); % Default element (Carbon)
        else
            elementList = H_hr.elementL; % Use defined elements
        end

        Atom_name_ = repmat(elementList, [NRPTS_, 1]); % Expand element names

        % Position replication
        for i = 1:NRPTS_
            sites_ = [sites_; positions + vectorList(i, :)];
        end

        % Plot expanded structure
        ax = POSCAR_plot(Rm_, sites_, Atom_name_,...
            'vectorL', vectorList,...
            'ax', ax,...
            'scale', options.scale,...
            'atomscale', options.atomscale,...
            'TwoD', options.TwoD,...
            'fast', options.fast);

    case 'HOPPING' % Hopping connection plot
        

        % First plot NRPTS base
        ax = H_hr.show('NRPTS',...
            'scale', options.scale,...
            'ax', ax,...
            'atomscale', options.atomscale,...
            'TwoD', options.TwoD,'vectorL', options.vectorList);

        Rm_ = double(Rm_ * options.scale); % Scaled lattice
        
        H_hr = H_hr.rewrite(); % Change To List Type

        % Handle coefficient/numeric hopping differently
        if H_hr.coe
            if isempty(options.Select)
                % Plot all coefficients
                Plot_Hopping(H_hr.vectorL, H_hr.HcoeL, Rm_, H_hr.orbL,...
                    'ax', ax, 'TwoD', options.TwoD, 'cmap', options.cmap,...
                    'fast',options.fast);
            else
                % Filter selected coefficients
                SelectL = ismember(H_hr.HcoeL, options.Select);
                Plot_Hopping(H_hr.vectorL(SelectL, :), H_hr.HcoeL(SelectL),...
                    Rm_, H_hr.orbL, 'ax', ax, 'TwoD', options.TwoD,...
                    'cmap', options.cmap,'fast',options.fast);
            end
        else
            if isempty(options.Select)
                % Plot all numeric hoppings
                %SelectL = H_hr.HnumL == options.Select;
                Plot_Hopping(H_hr.vectorL, H_hr.HnumL, Rm_, H_hr.orbL,...
                    'ax', ax, 'TwoD', options.TwoD);
            else
                % Filter numeric selections
                switch length(options.Select)
                    case 1 
                        SelectL = abs(H_hr.HnumL) > options.Select;
                    case 2
                        SelectL = all(H_hr.vectorL(:,4:5) == options.Select,2);
                    case 3
                        SelectL = all(H_hr.vectorL(:,1:3) == options.Select,2);
                end
                
                Plot_Hopping(H_hr.vectorL(SelectL, :), H_hr.HnumL(SelectL),...
                    Rm_, H_hr.orbL, 'ax', ax, 'TwoD', options.TwoD,...
                    'fast',options.fast);
            end
        end

        % 2D view configuration
        if options.TwoD
            view(ax, 0, 90); % Set 2D projection
            l = light(ax); % Add lighting
            l.Position = [0, 0, 1]; % Light direction
        end

    otherwise
        % Invalid mode handling (covered by arguments validation)
end

% Final title setting
title(options.ax, options.Title);
end