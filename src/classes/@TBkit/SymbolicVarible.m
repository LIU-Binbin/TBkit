function Symble = SymbolicVarible(seeds,superscript,subscript,level)
%SYMBOLICVARIBLE Create symbolic variable with formatted name
%   SYMBLE = SYMBOLICVARIBLE(SEEDS, SUPERSCRIPT, SUBSCRIPT, LEVEL) creates
%   a symbolic variable with formatted name including sub/superscripts.
%
%   Inputs:
%       seeds      - Base name string
%       superscript - Numeric superscript values
%       subscript  - Numeric subscript values
%       level      - (optional) Additional level index
%
%   Output:
%       Symble     - Created symbolic variable
%
%   Note:
%       Formats negative indices with '_bar' suffix
%       Adds '_ubar' suffix if level is provided
arguments
    seeds 
    superscript 
    subscript 
    level = [];
end
Superscript = "";
Subscript = "";
for i = 1:length(superscript)
    if superscript(i) < 0
        Superscript = Superscript+"__"+num2str(round(abs(superscript(i))))+"_bar";
    else
        Superscript = Superscript+"__"+num2str(round(abs(superscript(i))));
    end
end
for i = 1:length(subscript)
    if subscript(i) < 0
        Subscript = Subscript+"_"+num2str(round(abs(subscript(i))))+"_bar";
    else
        Subscript = Subscript+"_"+num2str(round(abs(subscript(i))));
    end
end
if ~isempty(level)
    Level = "_"+num2str(level)+"_ubar";
    Symble = sym(seeds+Subscript+Superscript+Level,'real');
else
    Symble = sym(seeds+Subscript+Superscript,'real');
end
end