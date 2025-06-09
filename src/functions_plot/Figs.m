function [fig, axes] = Figs(m,n,options)
arguments
    m = 1
    n = 1
    options.FontName = 'Helvetica'
    options.FontSize
    options.FigSize {mustBeMember(options.FigSize, ["auto", "standard", "wide"])} = "auto";
    options.Position = [];
end
%%
fig = figure();
set(fig, 'Color', 'w');
if isempty(options.Position)
    Unit = 'centimeters';
else
    Unit = 'normalized';
    Position = options.Position;
    set(fig,'unit',Unit ,...
    'position',Position );
end
switch options.FigSize
    case "auto"

        FontSize = 28;
    case "standard"
        fig_size = [10.5 11];
        fig_bias = [5 5];
        ax_size = [7 7]; 
        ax_bias = [2.8 2.6];

        FontSize = 24;
    case "wide"
        fig_size = [21 11];
        fig_bias = [5 5];
        ax_size = [7*2 7]; 
        ax_bias = [4.5 2.6];

        FontSize = 24;
end
if exist('options.FontSize','var')
    FontSize = options.FontSize;
end
%%


if options.FigSize ~= "auto"
    set(fig,'unit',Unit ,...
        'position', [fig_bias fig_size(1)*n fig_size(2)*m]);

end

p = 0;

axes = gobjects(m,n);
for i = 1:m
    for j = 1:n
        p = p+1;
        axes(i,j) = subplot(m,n,p);
        hold(axes(i,j), "on")

        if options.FigSize ~= "auto"
            set(axes(i,j), 'unit', 'centimeters', 'Position',...
                [ax_bias(1) + fig_size(1)*(j-1),...
                 ax_bias(2) + fig_size(2)*(m-i),...
                 ax_size]);
        else
            % axis(axes(i,j), 'square')
        end

        set(axes(i,j), 'LineWidth',1, ...
            'Box','on', ...
            'FontName', options.FontName, ...
            'FontSize', FontSize, ...
            'Layer', 'top')
    end
end
end