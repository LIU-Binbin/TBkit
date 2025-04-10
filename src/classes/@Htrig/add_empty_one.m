function H_htrig = add_empty_one(H_htrig,vector)
% ADD_EMPTY_ONE Adds empty terms to Htrig object with specified vector patterns.
%   H_htrig = ADD_EMPTY_ONE(H_htrig, vector) appends zero-valued terms to the 
%   Htrig object's coefficient/numerical lists or symbolic trigonometric list 
%   based on the input vector patterns. Skips duplicates for 'mat'/'list' types.
%
%   Inputs:
%       H_htrig - Target Htrig object to be modified.
%       vector  - Row vector(s) defining the hopping pattern. For 'mat'/'list' 
%                 types: each row represents a unique term. For other types: 
%                 symbolic vector defining the trigonometric component.
%
%   Outputs:
%       H_htrig - Updated Htrig object with new zero terms added.
%
%   Behavior by Type:
%       - 'mat'/'list': Adds zeros to HnumL/HcoeL matrices and registers vector 
%         in HsymL_numL/HsymL_coeL if not existing.
%       - Other types: Directly appends vector to HsymL_trig and initializes 
%         corresponding zero blocks in HcoeL/HnumL.
%
%   Example:
%       % Add empty term with vector [1, -1] to a trigonometric-type Htrig
%       H_htrig = add_empty_one(H_htrig, sym([1, -1]));
%
%   See also: Htrig, NumOrCoe

% Handle empty coefficient/numerical properties
if isempty(H_htrig.coe)||isempty(H_htrig.num)
    [~,~,H_htrig] = H_htrig.NumOrCoe();
end

switch H_htrig.Type
    case {'mat','list'}  % Matrix/List type handling
        % Numerical term processing
        if H_htrig.num
            vector = double(vector);
            for i = 1:size(vector,1)
                vector_single = (vector(i,:));
                try
                    % Skip existing terms
                    if (ismember(vector_single,H_htrig.HsymL_numL,'rows'))
                        continue;
                    end
                catch
                end
                % Add new term with zero matrix
                Kind = H_htrig.Kinds +1;
                H_htrig.HsymL_numL(Kind,:) = (vector_single);
                H_htrig.HnumL(:,:,Kind) = zeros(H_htrig.Basis_num);
            end
        end
        
        % Symbolic coefficient processing
        if H_htrig.coe
            vector = sym(vector);
            for i = 1:size(vector,1)
                vector_single = (vector(i,:));
                try
                    % Skip existing terms
                    if (ismember(vector_single,H_htrig.HsymL_coeL,'rows'))
                        continue;
                    end
                catch
                end
                % Add new term with symbolic zero matrix
                Kind = H_htrig.Kinds +1;
                H_htrig.HsymL_coeL(Kind,:) = (vector_single);
                H_htrig.HcoeL(:,:,Kind) = zeros(H_htrig.Basis_num,'sym');
            end
        end
        
    otherwise  % Trigonometric/other types
        Kind = H_htrig.Kinds+1;
        BASIS_NUM = H_htrig.Basis_num;
        % Directly append trigonometric term
        H_htrig.HsymL_trig(Kind) = vector;
        % Initialize corresponding zero blocks
        H_htrig.HcoeL(:,:,Kind) = zeros(BASIS_NUM,BASIS_NUM,'sym');
        H_htrig.HnumL(:,:,Kind)  = (zeros(BASIS_NUM,BASIS_NUM));
end
end