% DIFF_KLIST Differentiates a Htrig object along specified k-space directions.
%
% SYNTAX:
%   [varargout{1:nargout}] = diff_klist(H_htrig, dir, klist, options)
%
% DESCRIPTION:
%   This function differentiates a Htrig object with respect to one or more 
%   directions provided in the vector 'dir' (applicable for Htrig objects of type 'list').
%   For each specified direction, the function computes the derivative using the 
%   diff method and then converts the result to the HCAR representation using HCAR_gen,
%   with the k-points provided in 'klist'. Multiple outputs are returned, one for each
%   direction.
%
% INPUTS:
%   H_htrig   - A Htrig object.
%   dir       - A vector of integers specifying the differentiation directions (default: 1).
%   klist     - A list of k-points to be used in HCAR_gen conversion (default: []).
%   options   - A structure with field:
%                  Accuracy - A tolerance parameter for processing (default: 1e-6).
%
% OUTPUTS:
%   varargout - A cell array of Htrig objects in HCAR representation, corresponding
%               to the derivative along each direction in 'dir'. The number of outputs
%               is equal to the number of elements in 'dir'.
%
% EXAMPLE:
%   % Differentiate the Htrig object H along directions 1 and 2:
%   [H1, H2] = diff_klist(H, [1 2], klist, struct('Accuracy', 1e-6));
%
% SEE ALSO:
%   diff, HCAR_gen, namedargs2cell
%
function varargout = diff_klist(H_htrig, dir, klist, options)
arguments
    H_htrig
    dir = 1;
    klist = [];
    options.Accuracy = 1e-6;
end
if ~strcmp(H_htrig.Type, 'list')
end
optionscell = namedargs2cell(options);
H_htrig_tmp = H_htrig;
varargout{nargout} = H_htrig_tmp.diff(dir(numel(dir)), optionscell{:}).HCAR_gen(klist);
for i = 1:numel(dir)-1
    H_htrig_tmp = H_htrig;
    varargout{i} = H_htrig_tmp.diff(dir(i), optionscell{:}).HCAR_gen(klist);
end
end
