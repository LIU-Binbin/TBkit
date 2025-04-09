function str_out = quadk2K(str_in)
%QUADK2K Convert quadratic k-terms to compact notation
%
% Syntax:
%   str_out = quadk2K(str_in)
%
% Input:
%   str_in - Input string with k-terms (e.g., 'x^2')
%
% Output:
%   str_out - Converted string with compact notation (e.g., 'X')
%
% Description:
%   Standardizes k-space term notation by converting:
%   - x^2 → X, y^2 → Y, z^2 → Z
%   - Higher powers to combined forms (x^3 → X*x, etc.)
%
% Note:
%   Used for consistent symbolic representation
%   Preserves non-quadratic terms unchanged
%
% Example:
%   str = quadk2K('x^2 + y^3') % Returns 'X + Y*y'
str_out = strrep(str_in,'x^2',"X");
str_out = strrep(str_out,'y^2',"Y");
str_out = strrep(str_out,'z^2',"Z");
str_out = strrep(str_out,'x^3',"X*x");
str_out = strrep(str_out,'y^3',"Y*y");
str_out = strrep(str_out,'z^3',"Z*z");
str_out = strrep(str_out,'x^4',"X^2");
str_out = strrep(str_out,'y^4',"Y^2");
str_out = strrep(str_out,'z^4',"Z^2");
str_out = strrep(str_out,'x^5',"X^2*x");
str_out = strrep(str_out,'y^5',"Y^2*y");
str_out = strrep(str_out,'z^5',"Z^2*z");
end
