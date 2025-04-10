function dispAll(BasisFunction)
%DISPALL  Display all basis functions in a formatted manner.
%
%   DISPALL(BasisFunction) prints the details of each BasisFunc object in the 
%   array BasisFunction to the command window. For each element, it displays its 
%   coefficient (coe), basis function list (BFuncL), and spin information (spin).
%
%   Input:
%       BasisFunction - An array of BasisFunc objects.
%
%   Example:
%       dispAll(BasisFunction);
%
%   See also: BasisFunc

for i = 1:size(BasisFunction, 1)
    fprintf('============================================\n');
    % Display the first element in the current row.
    fprintf('%s', string(BasisFunction(i, 1).coe));
    fprintf('%s', string(BasisFunction(i, 1).BFuncL));
    fprintf('%s', string(BasisFunction(i, 1).spin));
    % Loop through remaining elements in the row.
    for j = 2:size(BasisFunction, 2)
        fprintf(' + ');
        fprintf('%s', string(BasisFunction(i, j).coe));
        fprintf('%s', string(BasisFunction(i, j).BFuncL));
        fprintf('%s', string(BasisFunction(i, j).spin));
    end
    fprintf('\n');
end
fprintf('============================================\n');
end

