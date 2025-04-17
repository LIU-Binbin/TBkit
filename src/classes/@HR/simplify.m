function H_hr = simplify(H_hr, Accuracy)
%SIMPLIFY Optimize HR object by removing insignificant terms
%   This function streamlines the HR object by:
%   - Removing terms below specified accuracy threshold
%   - Simplifying symbolic expressions
%   - Reducing vector basis dimensions through row reduction
%   - Maintaining data format consistency during cleanup
%
%   Inputs:
%       H_hr     : HR object containing Hamiltonian data
%       Accuracy : Threshold for term removal (default: 1e-6)
%                  Terms with absolute value below this threshold are removed
%
%   Output:
%       H_hr     : Simplified HR object with optimized storage
%
%   Features:
%       - Handles both numeric and symbolic coefficient types
%       - Supports list/matrix storage formats
%       - Processes vector hopping terms using row reduction
%       - Preserves structure while removing negligible elements

    %% Input validation and initialization
    if nargin < 2
        Accuracy = 1e-6;  % Set default precision threshold
    end
    
    %% Vector hopping specialization
    if H_hr.vectorhopping
        process_vector_basis();
        return;  % Early return for vector hopping case
    end
    
    %% Main simplification logic
    if H_hr.coe
        handle_symbolic_coefficients();
    end
    
    if H_hr.num
        handle_numeric_coefficients();
    end

    %% Nested functions for modular processing
    function process_vector_basis()
        % Simplify vector basis using row reduction and thresholding
        [H_hr.AvectorL, H_hr.BvectorL, H_hr.CvectorL] = deal(...
            process_single_basis(H_hr.AvectorL),...
            process_single_basis(H_hr.BvectorL),...
            process_single_basis(H_hr.CvectorL));
        
        function cleanedBasis = process_single_basis(basisMatrix)
            % Row-reduce and threshold individual basis matrix
            reducedBasis = rref(basisMatrix.').';  % Transpose for column operations
            rankEstimate = rank(reducedBasis);     % Get numerical rank
            cleanedBasis = real(reducedBasis(:,1:rankEstimate));  % Keep significant columns
            cleanedBasis(abs(cleanedBasis) < Accuracy) = 0;       % Threshold small values
        end
    end

    function handle_symbolic_coefficients()
        % Simplify symbolic expressions and remove zeros
        H_hr.HcoeL = simplify(H_hr.HcoeL);  % Algebraic simplification
        filter_zero_terms(H_hr.HcoeL);       % Apply thresholding
        
        function filter_zero_terms(coefficients)
            % Remove symbolic terms below threshold
            % zeroThreshold = sym(Accuracy);
            if strcmp(H_hr.Type, 'list')
                significantTerms = coefficients ~= sym(0);
                H_hr = H_hr.reseq(':', significantTerms);
            end
        end
    end

    function handle_numeric_coefficients()
        % Clean numerical coefficients based on accuracy
        switch H_hr.Type
            case 'list'
                filter_list_format();
            case 'mat'
                filter_matrix_format();
        end
        
        function filter_list_format()
            % Direct thresholding for list storage
            significantIndices = abs(H_hr.HnumL) > Accuracy;
            H_hr = H_hr.reseq(':', significantIndices);
        end
        
        function filter_matrix_format()
            % Matrix-wise thresholding for 3D storage
            toleranceMatrix = Accuracy * ones(size(H_hr.HnumL(:,:,1)));
            significantSlices = false(H_hr.NRPTS, 1);
            
            % Vectorized slice processing
            for sliceIdx = 1:H_hr.NRPTS
                currentSlice = H_hr.HnumL(:,:,sliceIdx);
                significantSlices(sliceIdx) = any(abs(currentSlice) > toleranceMatrix, 'all');
            end
            
            H_hr = H_hr.reseq(':', significantSlices);
        end
    end
end