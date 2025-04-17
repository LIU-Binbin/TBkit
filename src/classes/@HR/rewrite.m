function H_hr = rewrite(H_hr, options)
%REWRITE Convert HR object between matrix and list storage formats
%   This function reorganizes the Hamiltonian data storage format:
%   - Converts from 3D matrix format to list format when Type is not 'list'
%   - Converts back to matrix format when Type is 'list' and rewind option is true
%
%   Inputs:
%       H_hr    : HR object containing Hamiltonian data
%       options : Struct with processing parameters:
%           - rewind   : Flag to reverse conversion (default: false)
%           - Accuracy : Numerical precision for simplification (default: 1e-6)
%           - type     : Reserved for future type specifications
%
%   Output:
%       H_hr    : Modified HR object with converted data format

    % Argument validation and default options setup
    arguments
        H_hr;
        options.rewind = false;
        options.Accuracy = 1e-6;
        options.type = '';
    end
    
    WANNUM = H_hr.WAN_NUM;
    
    % Conversion from matrix to list format
    if ~strcmp(H_hr.Type, 'list')
        handle_matrix_conversion();
    % Conversion from list to matrix format (when rewind enabled)    
    elseif strcmp(H_hr.Type, 'list') && options.rewind
        handle_list_conversion();
    end
    
    H_hr = H_hr.simplify(options.Accuracy);

    %% Nested functions for modular organization
    function handle_matrix_conversion()
        % Handle conversion from 3D matrix to list format
        if H_hr.num
            SizeHopping = size(H_hr.HnumL);
            process_numeric_data(SizeHopping);
        elseif H_hr.coe && ~H_hr.num
            SizeHopping = size(H_hr.HcoeL);
            process_symbolic_data(SizeHopping);
        else
            initialize_empty_data();
        end
        H_hr.Type = 'list';
    end

    function handle_list_conversion()
        % Handle conversion from list format back to 3D matrix
        [vectorList, ~, icL] = unique(H_hr.vectorL(:,1:H_hr.Dim), 'rows');
        NRPTS_ = size(vectorList, 1);
        
        if H_hr.num
            HnumLtmp = zeros(WANNUM, WANNUM, NRPTS_);
            process_numeric_rewind(HnumLtmp, icL);
        else
            HcoeLtmp = sym(zeros(WANNUM, WANNUM, NRPTS_));
            process_symbolic_rewind(HcoeLtmp, icL);
        end
        
        H_hr.Type = 'mat';
        H_hr.vectorL = vectorList;
    end

    %% Data processing subfunctions
    function process_numeric_data(SizeHopping)
        if isvector(H_hr.HnumL)
            warning('May not need to rewrite. Do nothing');
            return;
        end
        % Reshape 3D numeric array to column vector
        H_hr.HnumL = reshape(H_hr.HnumL, [], 1);
        generate_vector_list(SizeHopping);
        
        if H_hr.overlap
            process_overlap_data(@H_hr.SnumL, 'vectorL_overlap');
        end
    end

    function process_symbolic_data(SizeHopping)
        if isvector(H_hr.HcoeL)
            warning('May not need to rewrite. Do nothing');
            return;
        end
        % Reshape 3D symbolic array to column vector
        H_hr.HcoeL = reshape(H_hr.HcoeL, [], 1);
        generate_vector_list(SizeHopping);
        
        if H_hr.overlap
            process_overlap_data(@H_hr.ScoeL, 'vectorL_overlap');
        end
    end

    %% Helper functions
    function generate_vector_list(SizeHopping)
        % Generate vector list from 3D matrix indices
        [iL, jL, nL] = ind2sub(SizeHopping, 1:prod(SizeHopping));
        H_hr.vectorL = [H_hr.vectorL(nL,:), iL.', jL.'];
    end

    function process_overlap_data(dataField, vectorField)
        % Process overlap data using function handles
        overlapData = dataField();
        sizeCoeL = size(overlapData);
        dataField(reshape(overlapData, [], 1));
        
        [iL, jL, nL] = ind2sub(sizeCoeL, 1:prod(sizeCoeL));
        H_hr.(vectorField) = [H_hr.(vectorField)(nL,:), iL.', jL.'];
    end

    function process_numeric_rewind(targetArray, indexConverter)
        % Convert numeric list back to 3D matrix
        iL = double(H_hr.vectorL(:,H_hr.Dim+1));
        jL = double(H_hr.vectorL(:,H_hr.Dim+2));
        indL = sub2ind(size(targetArray), iL, jL, indexConverter);
        targetArray(indL) = H_hr.HnumL;
        H_hr.HnumL = targetArray;
    end

    function process_symbolic_rewind(targetArray, indexConverter)
        % Convert symbolic list back to 3D matrix
        iL = double(H_hr.vectorL(:,H_hr.Dim+1));
        jL = double(H_hr.vectorL(:,H_hr.Dim+2));
        indL = sub2ind(size(targetArray), iL, jL, indexConverter);
        targetArray(indL) = H_hr.HcoeL;
        H_hr.HcoeL = targetArray;
    end

    function initialize_empty_data()
        H_hr.HnumL = [];
        H_hr.HcoeL = sym([]);
        H_hr.vectorL = [];
    end
end