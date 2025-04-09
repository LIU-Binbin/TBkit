function [fig, axes] = Figs(m,n,options)
arguments
    m = 1
    n = 1
    options.FontName = 'Helvetica'
    options.FigSize {mustBeMember(options.FigSize, ["large", "normal"])} = "large"
end
%%
fig_size = [10.5 9.9]; % 3x2 A4 size
fig_bias = [5 5];
ax_size = [7 7]; 
ax_bias = [2.8 1.8];

if options.FigSize == "large"
    fig_size = fig_size*2; % 3x2 A4 size
    fig_bias = fig_bias*2;
    ax_size = ax_size*2;
    ax_bias = ax_bias*2;
elseif options.FigSize == "normal"

end
%%
fig = figure();
set(fig,'unit', 'centimeters',...
    'position', [fig_bias fig_size(1)*n fig_size(2)*m],...
    'Color', 'w');

p = 0;
% alphabet_list = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n","o"];

axes = gobjects(m,n);
for i = 1:m
    for j = 1:n
        p = p+1;
        axes(i,j) = subplot(m,n,p);
        hold(axes(i,j), "on")
        set(axes(i,j), 'unit', 'centimeters', 'Position',...
            [ax_bias(1) + fig_size(1)*(j-1),...
             ax_bias(2) + fig_size(2)*(m-i),...
             ax_size]);
        set(axes(i,j), 'LineWidth',1,'Box','on', 'FontName', options.FontName)

        if options.FigSize == "large"
            set(axes(i,j), 'FontSize', 24);
        elseif options.FigSize == "normal"
            set(axes(i,j), 'FontSize', 18);
        end
        % title(axes(i,j), "("+alphabet_list(p)+")", 'Position',[-0.32, 1.05], 'Units','centimeters')
    end
end
end