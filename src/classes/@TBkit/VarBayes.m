function Varbayes = VarBayes(Varlist,VarGuess,VarWidth,VarFix)
% VARBAYES Create optimizable variables for Bayesian optimization
%
% Syntax:
%   Varbayes = VarBayes(Varlist,VarGuess,VarWidth)
%   Varbayes = VarBayes(Varlist,VarGuess,VarWidth,VarFix)
%
% Description:
%   Generates optimizableVariable objects for use with MATLAB's Bayesian
%   optimization routines.
%
% Input Arguments:
%   Varlist - List of variable names
%   VarGuess - Initial guesses for variables
%   VarWidth - Search width around initial guesses
%   VarFix - Variables to fix (optional)
%
% Output Arguments:
%   Varbayes - Array of optimizableVariable objects
%
% Example:
%   vars = VarBayes(["a","b"],[1,2],[0.5,0.5]);
if nargin < 4
    VarFix = sym([]);
end
for i = 1:length(Varlist)
    if isnan(VarGuess(i))
        VarGuess(i) = 0;
    end
    if ismember(Varlist(i),VarFix)
        Varbayes(i) = optimizableVariable(string(Varlist(i)),...
            [-VarWidth(i)+VarGuess(i),VarGuess(i)+VarWidth(i)],...
            'Optimize',true);
    else
        Varbayes(i) = optimizableVariable(string(Varlist(i)),...
            [-VarWidth(i)+ VarGuess(i),VarGuess(i)+VarWidth(i)]);
    end
end
end