classdef HollowKnight %(Abstract) A base class for hollow objects, extending its functionality in subclasses.
    % This class serves as a base class for objects that can be considered "hollow".
    % Subclasses of HollowKnight should define properties 'coe' and 'hollow'.

    properties
        coe {mustBeA(coe,["double","sym"])} = 1  % Coefficient, which can be a number or symbolic
    end

    properties(Abstract, Hidden)
        hollow; % Abstract property for the hollow state of the object
    end

    methods
        function obj = HollowKnight()
            % Constructor for HollowKnight, initializes the object.
            % This is an empty constructor to allow subclass initialization.
        end
    end
    
    %% Overload Methods for Mathematical Operations
    methods
        function C = mtimes(A,B)
            % only care structure
            if isa(A,'HollowKnight') && isa(B,'HollowKnight')
                if length(A) == 1 && length(B) == 1
                    C = A.*B;
                    return
                end
                % Define mtimes is very difficult. hard coding
                if isrow(A) && isrow(B)
                    C = A(1,1).*B(1,1);
                    for j = 2:size(B,2)
                        D = A(1).*B(j);
                        if ~D.hollow
                            C = [C,D];
                        end
                    end
                    for i = 2:size(A,2)
                        for j = 1:size(B,2)
                            D = A(i).*B(j);
                            if ~D.hollow
                                C = [C,D];
                            end
                        end
                    end
                    C = contractrow(C);
                    return;
                end
                % 1*col = col multipliy means
                if isrow(A) && ~isrow(B)
                    C = A*B(1,:);
                    for  j= 2:size(B,1)
                        C = [C;A*B(j,:)];
                    end
                    return;
                end
                % col*col  = col*col multipliy means kron
                if ~isrow(A) && ~isrow(B)
                    C = A(1,:)*B;
                    for i = 2:size(A,1)
                        C = [C;A(i,:)*B];
                    end
                    return;
                end
            elseif ~isa(A,'HollowKnight') && isa(B,'HollowKnight')
                % Define mtimes is very difficult. hard coding
                if length(A) == 1 && length(B) == 1
                    C = A.*B;
                    return
                end
                % 1 *(m*n)
                if isscalar(A)
                    C = B;
                    for i = 1:numel(C)
                        C(i) = A.*B(i);
                    end
                    C = contract(C);
                    return
                end
                % row * row
                if isrow(A) && isrow(B)
                    if length(A) == length(B)
                        C = A(1).*B(1);
                        for j = 2:size(B,2)
                            D = A(j).*B(j);
                            if ~D.hollow
                                C = [C,D];
                            end
                        end
                        C = contractrow(C);
                        return;
                    else
                        error('not imply yet')
                    end
                end
                % row*col = 1*? multipliy means
                if isrow(A) && ~isrow(B)
                    if size(A,2) == size(B,1)
                        C = A(1)*B(1,:);
                        for  j= 2:size(B,1)
                            C = [C,A(j)*B(j,:)];
                        end
                        C = contractrow(C);
                    else
                        error('not imply yet')
                    end
                    return;
                end
                % (row*col)*col  = row*? multipliy means kron
                if ~isrow(A) && ~isrow(B)
                    if size(A,2) == size(B,1)
                        C = A(1,:)*B;
                        for i = 2:size(A,1)
                            C = [C;A(i,:)*B];
                        end
                    else
                        error('not imply yet')
                    end
                    return;
                end
            elseif isa(A,'HollowKnight') && ~isa(B,'HollowKnight')
                C = mtimes(B,A);
            else
            end
        end
        
        function C = times(A,B)
            if isa(A,'HollowKnight') && isa(B,'HollowKnight')
                % return times
                C = innertimes(A,B);
            elseif ~isa(A,'HollowKnight') && isa(B,'HollowKnight')
                C = B;
                if length(A) == 1
                    for i =1:numel(B)
                        C(i).coe = A.*C(i).coe ;
                    end
                elseif isvector(A) && length(A) == size(B,1)
                    for i = 1:size(B,1)
                        C(i,:) = A(i).*B(i,:);
                    end
                elseif isequal(size(A) ,size(B))
                    for i = 1:size(B,1)
                        for j = 1:size(B,2)
                            C(i,j) = A(i,j).*B(i,j);
                        end
                    end
                else
                    error('times wrong')
                end
            elseif isa(A,'HollowKnight') && ~isa(B,'HollowKnight')
                C = times(B,A);% suppose commute!!
            else

            end
        end

        function C = innertimes(A, B)
            % Inner product of two HollowKnight objects.
            % Computes the inner product of two HollowKnight objects, multiplying their coefficients.
            C = A;
            C.coe = A.coe * B.coe;
        end
    end

    %% Rotation Methods
    methods
        function Am = rotate(A, rotm, rightorleft, options)
            % Rotates HollowKnight objects based on a rotation matrix or quaternion.
            % Rotates the object using a specified rotation matrix or quaternion.
            
            arguments
                A HollowKnight;         % HollowKnight object
                rotm {mustBeSize(rotm, [3, 3; 1, 5;1 4])} = diag([1, 1, 1]); % Rotation matrix or quaternion
                rightorleft = 'right'; % Rotation direction
                options.sym = true;     % Whether to use symbolic rotation
                options.conjugate = false; % Whether to use conjugate symmetry
                options.antisymmetry = false; % Antisymmetry flag
            end
            
            optionsCell = namedargs2cell(options);
            Am = rotaterow(A(1,:), rotm, rightorleft, optionsCell{:});
            for i = 2:size(A, 1)
                Am = [Am; rotaterow(A(i,:), rotm, rightorleft, optionsCell{:})];
            end
        end
        
        function A_Lj = rotaterow(A, rotm, rightorleft, options)
            % Rotate each row of the HollowKnight object.
            % Rotates each row in the given HollowKnight object using the specified rotation matrix.
            
            arguments
                A HollowKnight; 
                rotm {mustBeSize(rotm, [3, 3; 1, 5;1 4])} = diag([1, 1, 1]); 
                rightorleft = 'right';
                options.sym = true;
                options.conjugate = false;
                options.antisymmetry = false;
            end
            
            optionsCell = namedargs2cell(options);
            A_Lj = contractrow(rotatesingle(A(1), rotm, rightorleft, optionsCell{:}));
            for i = 2:size(A, 2)
                A_Lj = [A_Lj, contractrow(rotatesingle(A(i), rotm, rightorleft, optionsCell{:}))];
            end
        end
         function Ak = rotatesingle(A,rotm,rightorleft,options)
            arguments
                A HollowKnight;
                rotm {mustBeSize(rotm,[3 3;1 5;5 1;1 4;4 1])}= diag([1 1 1]);% [alpha beta gamma]
                rightorleft = 'right';
                options.sym = true;
                options.conjugate = false;
                options.antisymmetry = false;
            end
            optionsCell = namedargs2cell(options);
            if strcmp(rightorleft,'right')
                RightorLeft = 1; % Check which is right!!!!!
            else
                RightorLeft = -1;
            end
            if isequal(size(rotm),[3 3])
                if det(rotm) == -1
                    % inversion
                    immproper = true;
                    rotm = -rotm;
                else
                    immproper = false;
                end
                %
                % note the robot toolbos use left axis! clockwise rotation
                % here we change our right rotm to left rotm
                rotm = inv(rotm);
                abc = Oper.Rotation2eul(rotm);% alpha beta gamma in ZYZ
            elseif isequal(size(rotm),[4 1]) || isequal(size(rotm),[1 4])
                if sym(abs(rotm(end))) ~= sym(1)
                    warning('wrong input,eular angle ZYZ right format:[alpha beta gamma det()]');
                    abc = Oper.axang2eul(rotm(1:4));
                    immproper = false;
                else
                    if sym(rotm(end)) == sym(-1)
                        immproper = true;
                    else
                        immproper = false;
                    end
                    abc = rotm(1:3);
                end
            elseif isequal(size(rotm),[5 1]) || isequal(size(rotm),[1 5])
                if sym(abs(rotm(end))) ~= sym(1)
                    error('wrong input,axis angle det right format:[nx ny nz theta det()]');
                end
                if sym(rotm(end)) == sym(-1)
                    immproper = true;
                else
                    immproper = false;
                end
                abc = Oper.axang2eul(rotm(1:4));
            end
            if options.sym
                abc = sym(abc);
            else
                %abc
            end
            Ak = rotateinner(A,abc,RightorLeft,immproper,options.conjugate,options.antisymmetry);
        end
    
    end

    %% Contract Methods
    methods
        function A = horzcat(varargin)
            % Horizontal concatenation of HollowKnight objects.
            % Concatenate multiple HollowKnight objects along the horizontal axis.
            
            arguments(Repeating)
                varargin HollowKnight;
            end
            
            for i = 1:nargin
                classL(i) = string(class(varargin{i}));
            end
            classname = string(class(varargin{1}));
            chooseL = find(classname ~= classL);
            
            for i = 1:chooseL
                varargin{i} = feval(classname, varargin{i});
            end
            A = builtin('horzcat', varargin{:});
        end
        
        function A = vertcat(varargin)
            % Vertical concatenation of HollowKnight objects.
            % Concatenate multiple HollowKnight objects along the vertical axis.
            
            arguments(Repeating)
                varargin HollowKnight;
            end
            
            for i = 1:nargin
                ncol(i) = size(varargin{i}, 2);
            end
            ncolmax = max(ncol);
            for i = 1:nargin
                if ncolmax ~= ncol(i)
                    varargin{i} = varargin{i}.padded(ncolmax);
                end
            end
            A = builtin('vertcat', varargin{:});
        end
        
        function A = padded(A, ncol, nrow)
            % Pad the HollowKnight object to match the desired number of columns.
            % Adds columns or rows if needed.
            
            arguments
                A HollowKnight;
                ncol = 1;
                nrow = 1;
            end
            
            if ncol - size(A, 2) <= 0
                return;
            end
            nrow = size(A, 1);
            A = [A, HollowKnight.creator(ncol - size(A, 2), nrow, class(A))];
        end
        
        function A = cleanrow(A)
            % Clean up a row by removing empty (hollow) elements.
            % This removes elements with NaN coefficients from the row.
            
            if isrow(A)
                coeList = [A.coe];
                coeList(isnan(coeList)) = 0;
                A(sym(coeList) == sym(0)) = [];
            else
                error('Not a row!');
            end
        end
        
        function B = contract(A, options)
            % Contract a HollowKnight object along its rows.
            % This method reduces dimensions by merging rows based on certain options.
            
            arguments
                A HollowKnight;
                options.forgetcoe = false;
            end
            
            optionsCell = namedargs2cell(options);
            nrow = size(A, 1);
            B = contractrow(A(1, :));
            for i = 2:nrow
                HollowKnightRow = contractrow(cleanrow(A(i, :)), optionsCell{:});
                B = [B; HollowKnightRow];
            end
        end
    end

    methods(Abstract)
        B = contractrow(A);
    end
    
    methods
        function A = HollowMe(A)
            % Marks the object as a "hollow" object by setting the coefficient to NaN.
            % This method makes the object hollow by setting its coefficient to NaN.
            try
                A.coe = nan;
            catch
                for i = numel(A)
                    A(i).coe = nan;
                end
            end
            
        end
    end

    %% Static Helper Methods
    methods(Static)
        function [comparerow_unique,sumrow_unique] = generalcontractrow(comparerow,sumrow)
            arguments
                comparerow
                sumrow
            end
            % rm nan
            selecL =  ~logical(sum(isnan(sumrow),2));
            %
            [comparerow_unique,sumrow_unique] = HollowKnight.generalcontractrow2(comparerow(selecL,:),sumrow(selecL,:));
            Lic = find(sum(logical(zeros(1,size(sumrow,2),class(sumrow))==sumrow_unique),2));
            sumrow_unique(Lic,:) = [];
            comparerow_unique(Lic,:) = [];
        end

        function [comparerow_unique, sumrow_unique] = generalcontractrow2(comparerow, sumrow)
            % This function contracts rows based on unique comparisons and sums them.
            % It returns unique comparerows and the corresponding sumrow_unique.
            %
            % Arguments:
            % - comparerow: Array of rows to compare.
            % - sumrow: Array of values to sum for each unique comparerow.
            %
            % Returns:
            % - comparerow_unique: The unique rows from comparerow.
            % - sumrow_unique: The sum of the rows corresponding to each unique comparerow.

            arguments
                comparerow
                sumrow
            end
            
            % Case when comparerow has only one row
            if size(comparerow, 1) == 1
                comparerow_unique = comparerow;
                sumrow_unique = sumrow;
                return;
            end
            
            % Get unique rows in comparerow
            [comparerow_unique, ~, ~] = unique(comparerow, 'rows');
            
            % If all rows are unique, return input as is
            if size(comparerow_unique, 1) == size(comparerow, 1)
                comparerow_unique = comparerow;
                sumrow_unique = sumrow;
                return;
            end
            
            % Initialize sumrow_unique with zeros
            sumrow_unique = zeros([size(comparerow_unique, 1), size(sumrow, 2)], class(sumrow));
            
            % Loop through unique comparerow and sum corresponding rows in sumrow
            for i = 1:size(comparerow_unique, 1)
                comparerow_unique_tmp = comparerow_unique(i, :);
                [~, index_ic] = ismember(comparerow, comparerow_unique_tmp, 'rows');
                sumrow_unique(i, :) = sum(sumrow(logical(index_ic), :), 1);
            end
        end
        function HollowKnightObj = creator(ncol, nrow, classname)
            % Creates a HollowKnight object and initializes it as hollow.
            % It fills a matrix of size nrow x ncol with the hollow object.
            %
            % Arguments:
            % - ncol: Number of columns to create.
            % - nrow: Number of rows to create.
            % - classname: The class name for creating the HollowKnight object.
            %
            % Returns:
            % - HollowKnightObj: A matrix of hollow objects of specified size.

            HollowObj = feval(classname);  % Create an instance of the specified class
            HollowObj = HollowObj.HollowMe();  % Make the object hollow
            HollowKnightObj = repmat(HollowObj, [nrow, ncol]);  % Create a matrix of hollow objects
        end
        function varargout = StandardInput(varargin, options)
            % This function automatically reshapes the input arguments into matrices
            % of the specified number of rows and columns.
            %
            % Arguments:
            % - varargin: Input arguments that may be scalars, vectors, or matrices.
            % - options: Structure containing 'nrow' and 'ncol' to specify the number of rows and columns.
            %
            % Returns:
            % - varargout: The reshaped inputs.

            arguments(Repeating)
                varargin
            end
            arguments
                options.nrow = 1;  % Default number of rows
                options.ncol = 1;  % Default number of columns
            end
            
            % Extract the number of rows and columns from options
            nrow = options.nrow;
            ncol = options.ncol;
            
            % Loop through each input argument and reshape
            for i = 1:numel(varargin)
                % If the input is a scalar, expand it to the matrix size
                if isscalar(varargin{i})
                    tmp = repmat(varargin{i}, [nrow, ncol]);
                else
                    tmp = varargin{i};
                end
                
                % If the input is a column vector, reshape to nrow x ncol
                if iscolumn(tmp)
                    tmp = repmat(tmp, [1, ncol]);
                % If the input is a row vector, reshape to nrow x ncol
                elseif isrow(tmp)
                    tmp = repmat(tmp, [nrow, 1]);
                end
                
                % Store the reshaped input into varargout
                varargout{i} = tmp;
            end
        end
    end
    %% operator
    methods
    % InnerProduct
    function Smat = InnerProduct(A, B, options)
        % Calculate the inner product matrix between two HollowKnight objects A and B.
        % Options can specify symmetry, strictness, union, and squareness.
        arguments
            A HollowKnight
            B HollowKnight
            options.sym = false;       % Whether to use symbolic computation
            options.strict = false;    % Whether to enforce strict matching
            options.union = false;     % Whether to perform union operation
            options.square = true;     % Whether to enforce square matrix
        end
        optionsCell = namedargs2cell(options);
        if size(A, 1) ~= size(B, 1)
            % Attempt to contract A and B if their sizes do not match
            A = contract(A);
            B = contract(B);
            if length(A) ~= length(B) && options.union
                % Perform union operation if specified
                % (Implementation needed)
            end
        end
        % Initialize the result matrix with appropriate data type
        if options.sym
            Ssingle = sym(0);
        else
            Ssingle = 0;
        end
        Smat = repmat(Ssingle, [size(A, 1), size(B, 1)]);
        for i = 1:size(A, 1)
            for j = 1:size(B, 1)
                % Compute the inner product for each pair of rows
                Smat(i, j) = InnerProduct_row(A(i, :), B(j, :), optionsCell{:});
            end
        end
    end
    function SingleSum = InnerProduct_row(A_row, B_row, options)
        % Calculate the inner product for a pair of rows from HollowKnight objects A and B.
        % Options can specify symmetry, strictness, union, and squareness.
        arguments
            A_row HollowKnight
            B_row HollowKnight
            options.sym = true;        % Whether to use symbolic computation
            options.strict = false;    % Whether to enforce strict matching
            options.union = false;     % Whether to perform union operation
            options.square = true;     % Whether to enforce square matrix
        end
        % Initialize the sum with appropriate data type
        if options.sym
            SingleSum = sym(0);
        else
            SingleSum = 0;
        end
        if isrow(A_row) && isrow(B_row)
            for i = 1:size(A_row, 2)
                iA_row = A_row(i);
                if iA_row.hollow
                    continue;
                end
                for j = 1:size(B_row, 2)
                    jB_row = B_row(j);
                    if jB_row.hollow
                        continue;
                    end
                    if iA_row == jB_row
                        % Accumulate the product of coefficients for matching elements
                        SingleSum = SingleSum + A_row(i).coe * B_row(j).coe;
                    end
                end
            end
        end
    end
    function HollowKnightObj = Prow(HollowKnightObj, P, options)
        % Permute the rows of the HollowKnight object based on permutation vector P.
        % Options can specify whether the permutation is antisymmetric.
        arguments
            HollowKnightObj
            P
            options.antisymmetric = true; % Whether the permutation is antisymmetric
        end
        if isvector(P)
            permVec = P;
        else
            permVec = permMatToVec(P);
        end
        % Determine the parity of the permutation
        if options.antisymmetric
            Parity = ParityPerm(permVec)^(HollowKnightObj(1).J * 2);
        else
            Parity = 1;
        end
        % Apply the permutation to the rows
        HollowKnightObj = Parity .* HollowKnightObj(:, permVec);
    end
    function HollowKnightObj = Pcol(HollowKnightObj, P, options)
        % Permute the columns of the HollowKnight object based on permutation vector P.
        % Options can specify whether the permutation is antisymmetric.
        arguments
            HollowKnightObj
            P
            options.antisymmetric = false; % Whether the permutation is antisymmetric
        end
        if isvector(P)
            permVec = P;
        else
            permVec = permMatToVec(P);
        end
        % Determine the parity of the permutation
        if options.antisymmetric
            Parity = ParityPerm(permVec)^(HollowKnightObj(1).J * 2);
        else
            Parity = 1;
        end
        % Apply the permutation to the columns
        HollowKnightObj = Parity .* HollowKnightObj(permVec, :);
    end
    function HollowKnightObj = Pmat(HollowKnightObj, P, options)
        % Permute the matrix of the HollowKnight object based on permutation matrix P.
        % Options can specify whether the permutation is antisymmetric.
        arguments
            HollowKnightObj HollowKnight
            P
            options.antisymmetric = false; % Whether the permutation is antisymmetric
        end
        if isvector(P)
            permMat = permVecToMat(P);
        else
            permMat = P;
        end
        % Apply the permutation to the matrix
        HollowKnightObj = permMat * HollowKnightObj;
    end
    end
 
end