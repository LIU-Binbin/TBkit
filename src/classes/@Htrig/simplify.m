function H_htrig = simplify(H_htrig, Accuracy, options)
%SIMPLIFY Simplify the hopping terms in a Htrig object.
%
%   H_htrig = SIMPLIFY(H_htrig)
%   performs symbolic simplification on the hopping matrices stored in a
%   Htrig object. It removes zero or near-zero terms based on a default
%   accuracy threshold (1e-6).
%
%   H_htrig = SIMPLIFY(H_htrig, Accuracy)
%   specifies the threshold for determining near-zero numerical values.
%
%   H_htrig = SIMPLIFY(H_htrig, Accuracy, options)
%   allows additional options. Available options:
%       options.reduce (default: false)
%           - If true, attempts to reduce variable terms to a canonical form
%             using TBkit.cleanVar based on the given Accuracy.
%
%   Inputs:
%       H_htrig  - An instance of the Htrig class
%       Accuracy - Numerical threshold below which values are treated as zero
%       options  - Struct with fields:
%                   reduce: logical, whether to apply symbolic variable reduction
%
%   Output:
%       H_htrig  - Updated Htrig object with simplified terms
%
%   Behavior:
%       - Both symbolic (`HcoeL`) and numerical (`HnumL`) terms are filtered.
%       - Representation-specific handling is performed based on Htrig.Type:
%         - 'list': sparse representation with vectorized storage
%         - 'mat', 'exp', 'sincos': full tensor storage
%       - Zero-valued slices are removed to reduce memory footprint and
%         improve performance.
%
%   See also: TBkit.cleanVar, simplify (Symbolic Math Toolbox), reseq

arguments
    H_htrig Htrig;
    Accuracy = 1e-6;
    options.reduce = false;
end

% Optional symbolic reduction
if options.reduce
    for i = 1:numel(H_htrig.HcoeL)
        if H_htrig.HcoeL(i) ~= sym(0)
            H_htrig.HcoeL(i) = TBkit.cleanVar(H_htrig.HcoeL(i), log10(Accuracy));
        end
    end
end

% Refresh numeric/symbolic status
[~, ~, H_htrig] = NumOrCoe(H_htrig);

% Simplify symbolic part
if H_htrig.coe
    H_coeL_tmp = simplify(H_htrig.HcoeL);
    H_htrig.HcoeL = H_coeL_tmp;
    switch H_htrig.Type
        case 'list'
            Kinds_list = find(H_coeL_tmp ~= sym(0));
            H_htrig = H_htrig.reseq(':', Kinds_list);
        case {'mat', 'exp', 'sincos'}
            zerosMat = zeros(size(H_coeL_tmp(:, :, 1)), 'sym');
            Kinds_list = true(H_htrig.Kinds, 1);
            for i = 1:H_htrig.Kinds
                if isequal(zerosMat, H_coeL_tmp(:, :, i))
                    Kinds_list(i) = false;
                end
            end
            H_htrig = H_htrig.reseq(':', Kinds_list);
    end
end

% Simplify numerical part
if H_htrig.num
    H_numL_tmp = H_htrig.HnumL;
    switch H_htrig.Type
        case 'list'
            Kinds_list = find(abs(H_numL_tmp) > Accuracy);
            H_htrig = H_htrig.reseq(':', Kinds_list);
        case {'mat', 'exp', 'sincos'}
            zerosMat = ones(H_htrig.Basis_num) * Accuracy;
            Kinds_list = true(H_htrig.Kinds, 1);
            for i = 1:H_htrig.Kinds
                if ~any(any(abs(H_numL_tmp(:, :, i)) > zerosMat))
                    Kinds_list(i) = false;
                end
            end
            H_htrig = H_htrig.reseq(':', Kinds_list);
    end
end
end

