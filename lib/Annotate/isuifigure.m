function tf = isuifigure(f)
%ISUIFIGURE Determine if handle references a uifigure.
%
% Usage:
%
%   TF = ISUIFIGURE(F)
%
% Inputs:
%
%   F <1x1 handle>
%     - Handle referencing a potential uifigure
%
% Outputs:
%
%   TF <1x1 logical>
%     - True if handle references a uifigure, false otherwise
%
% References:
%
%   https://www.mathworks.com/matlabcentral/answers/348387-distinguish-uifigure-from-figure-programmatically

    try
        tf = ~isempty(matlab.ui.internal.dialog.DialogHelper.getFigureID(f));
    catch
        try
            tf = matlab.ui.internal.isUIFigure(f);
        catch
            tf = false;
        end
    end
end
