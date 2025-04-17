function Line = Plot_line(Start, End, options)
%PLOT_LINE Plots a 3D line from Start to End with customizable properties.
%
%   Line = Plot_line(Start, End, options)
%
%   Inputs:
%       Start  - 1x3 vector [x, y, z], start point.
%       End    - 1x3 vector [x, y, z], end point.
%       options - structure with fields:
%           .ax           - Axes handle (default: gca)
%           .Color        - Line color (default: 'b')
%           .LineStyle    - Line style (default: '-')
%           .LineWidth    - Line width (default: 1.0)
%           .DisplayName  - Display name for legend (default: 'Hopping')
%           .W, .H        - Reserved fields (currently unused)
%           .CData        - Optional: used for coloring markers or other props
%
%   Output:
%       Line - Line object handle

    arguments
        Start (1,3) double = [0,0,0]
        End   (1,3) double = [1,1,1]
        options.ax = gca
        options.Color = 'b'
        options.LineStyle = '-'
        options.LineWidth (1,1) double = 1.0
        options.DisplayName = 'Hopping'
        options.W = []    % reserved
        options.H = []    % reserved
        options.CData = []  % reserved
    end

    % Draw the line using plot3
    Line = plot3(options.ax, ...
        [Start(1), End(1)], ...
        [Start(2), End(2)], ...
        [Start(3), End(3)], ...
        options.LineStyle, ...
        'Color', options.Color, ...
        'LineWidth', options.LineWidth, ...
        'DisplayName', options.DisplayName ...
    );
end
