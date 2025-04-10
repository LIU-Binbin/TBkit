% RESEQ Reorders the basis (Wannier functions) and/or hopping terms (kinds) of Htrig object.
%
% SYNTAX:
%   H_htrig = reseq(H_htrig, wan_list)
%   H_htrig = reseq(H_htrig, wan_list, kinds_list)
%
% DESCRIPTION:
%   This function allows the user to reorder or filter the Wannier orbitals and/or
%   hopping term "kinds" in an Htrig object, depending on the Hamiltonian storage type.
%
% INPUT:
%   H_htrig     - A Htrig object.
%   wan_list    - A list of indices specifying which Wannier orbitals to retain/reorder,
%                 or ':' to keep all.
%   kinds_list  - A list of indices specifying which hopping kinds to retain/reorder,
%                 or ':' to keep all. (Optional, defaults to ':')
%
% OUTPUT:
%   H_htrig     - The modified Htrig object with the specified reordering or filtering applied.
%
% BEHAVIOR BY TYPE:
%   - 'sparse':
%       HnumL is a cell array; elements are subset based on `wan_list` and `kinds_list`.
%
%   - 'list':
%       Data stored in vector and matrix list form. Indexing is more complex and uses `ismember`
%       logic for `wan_list`. Kinds are filtered directly via indexing.
%
%   - 'mat' or other types:
%       Assumes 3D matrix representation for both numeric and symbolic Hamiltonian terms.
%
% NOTES:
%   - If `wan_list` is not ':', all corresponding site/orbital/quantum-related fields are also updated.
%   - Empty symbolic arrays are replaced with `sym([])` for consistency.
%
% EXAMPLE:
%   % Keep only the first 4 Wannier orbitals and the first 10 kinds
%   H = reseq(H, 1:4, 1:10);
%
function H_htrig = reseq(H_htrig, wan_list, kinds_list)
    if nargin < 3
        kinds_list = ':';
    end

    % WAN_RESEQ (Wannier orbital dimension trimming)
    if ~isequal(wan_list, ':')
        switch H_htrig.Type
            case 'sparse'
                for i = 1:H_htrig.Kinds
                    H_htrig.HnumL{i} = H_htrig.HnumL{i}(wan_list, wan_list);
                end
            case 'list'
                mask = all(ismember(H_htrig.vectorL(:, H_htrig.Dim+1:H_htrig.Dim+2), int32(wan_list)), 2);
                if H_htrig.num
                    H_htrig.HnumL = H_htrig.HnumL(mask, :);
                    H_htrig.HsymL_numL = H_htrig.H_htrig.HsymL_numL(mask, :);
                    H_htrig.HsymL_coeL = [];
                    H_htrig.HcoeL = [];
                end
                if H_htrig.coe
                    H_htrig.HcoeL = H_htrig.HcoeL(mask, :);
                    H_htrig.HsymL_coeL = H_htrig.H_htrig.HsymL_coeL(mask, :);
                    H_htrig.HsymL_numL = [];
                    H_htrig.HnumL = [];
                end
            case 'mat'
                if H_htrig.num
                    H_htrig.HnumL = H_htrig.HnumL(wan_list, wan_list, :);
                else
                    H_htrig.HnumL = [];
                    H_htrig.HsymL_numL = [];
                end
                if H_htrig.coe
                    H_htrig.HcoeL = H_htrig.HcoeL(wan_list, wan_list, :);
                else
                    H_htrig.HcoeL = sym([]);
                end
            otherwise
                if H_htrig.num
                    H_htrig.HnumL = H_htrig.HnumL(wan_list, wan_list, :);
                else
                    H_htrig.HnumL = [];
                end
                if H_htrig.coe
                    H_htrig.HcoeL = H_htrig.HcoeL(wan_list, wan_list, :);
                else
                    H_htrig.HcoeL = sym([]);
                end
        end

        % Apply the same filtering to auxiliary data if available
        if ~isempty(H_htrig.sites)
            try
                H_htrig.sites = H_htrig.sites(wan_list);
            catch
            end
        end
        if ~isempty(H_htrig.orbL)
            H_htrig.orbL = H_htrig.orbL(wan_list, :);
        end
        if ~isempty(H_htrig.elementL)
            H_htrig.elementL = H_htrig.elementL(wan_list, :);
        end
        if ~isempty(H_htrig.quantumL)
            H_htrig.quantumL = H_htrig.quantumL(wan_list, :);
        end
    end

    % KIND_RESEQ (Hopping kinds trimming)
    if ~isequal(kinds_list, ':')
        switch H_htrig.Type
            case 'sparse'
                H_htrig.HnumL = H_htrig.HnumL(kinds_list);
            case 'list'
                if H_htrig.num
                    H_htrig.HsymL_numL = H_htrig.HsymL_numL(kinds_list, :);
                    H_htrig.HnumL = H_htrig.HnumL(kinds_list);
                    if ~isempty(H_htrig.HcoeL)
                        try
                            H_htrig.HcoeL = H_htrig.HcoeL(:, :, kinds_list);
                        catch
                        end
                    end
                    if ~isempty(H_htrig.HsymL_coeL)
                        try
                            H_htrig.HsymL_coeL = H_htrig.HsymL_coeL(kinds_list, :);
                        catch
                        end
                    end
                end
                if H_htrig.coe
                    H_htrig.HsymL_coeL = H_htrig.HsymL_coeL(kinds_list, :);
                    H_htrig.HcoeL = H_htrig.HcoeL(kinds_list);
                    if ~isempty(H_htrig.HnumL)
                        try
                            H_htrig.HnumL = H_htrig.HnumL(:, :, kinds_list);
                        catch
                        end
                    end
                    if ~isempty(H_htrig.HsymL_numL)
                        try
                            H_htrig.HsymL_numL = H_htrig.HsymL_numL(kinds_list, :);
                        catch
                        end
                    end
                end
            case 'mat'
                if H_htrig.num
                    H_htrig.HsymL_numL = H_htrig.HsymL_numL(kinds_list, :);
                    H_htrig.HnumL = H_htrig.HnumL(:, :, kinds_list);
                    if ~isempty(H_htrig.HcoeL)
                        try
                            H_htrig.HcoeL = H_htrig.HcoeL(:, :, kinds_list);
                        catch
                        end
                    end
                    if ~isempty(H_htrig.HsymL_coeL)
                        try
                            H_htrig.HsymL_coeL = H_htrig.HsymL_coeL(kinds_list, :);
                        catch
                        end
                    end
                end
                if H_htrig.coe
                    H_htrig.HsymL_coeL = H_htrig.HsymL_coeL(kinds_list, :);
                    H_htrig.HcoeL = H_htrig.HcoeL(:, :, kinds_list);
                    if ~isempty(H_htrig.HnumL)
                        try
                            H_htrig.HnumL = H_htrig.HnumL(:, :, kinds_list);
                        catch
                        end
                    end
                    if ~isempty(H_htrig.HsymL_numL)
                        try
                            H_htrig.HsymL_numL = H_htrig.HsymL_numL(kinds_list, :);
                        catch
                        end
                    end
                end
            otherwise
                H_htrig.HsymL_trig = H_htrig.HsymL_trig(kinds_list);
                if H_htrig.coe
                    H_htrig.HcoeL = H_htrig.HcoeL(:, :, kinds_list);
                end
                if H_htrig.num
                    H_htrig.HnumL = H_htrig.HnumL(:, :, kinds_list);
                end
        end
    end
end

