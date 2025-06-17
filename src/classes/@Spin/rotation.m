function U = rotation(SpinObj, rotm, rightorleft, options)
    % ROTATION applies a rotation to a spin object based on input rotation matrix or axis-angle format.
    % 
    % Input Arguments:
    %   SpinObj      - Spin object (of class Spin)
    %   rotm         - Rotation matrix (3x3), or axis-angle (4x1/1x4), or quaternion (5x1/1x5)
    %   rightorleft  - Rotation direction ('right' or 'left')
    %   options      - Optional named parameters:
    %     - full      - Whether to return a full matrix (default false)
    %     - sym       - Whether to use symbolic computation (default true)
    %
    % Output:
    %   U            - The resulting rotation operator (as a matrix)

    % Parsing optional arguments using MATLAB's 'arguments' block
    arguments
        SpinObj Spin;                                % Spin object (must be of class 'Spin')
        rotm {mustBeSize(rotm, [3 3; 1 5; 5 1; 1 4; 4 1])} = diag([1 1 1]); % Rotation matrix or axis-angle
        rightorleft = 'right';                       % Direction of rotation ('right' or 'left')
        options.full = false;                        % Option to return full matrix (default false)
        options.sym = true;                          % Use symbolic computation (default true)
    end
    
    % Convert named options to cell array for later use
    optionsCell = namedargs2cell(options);
    
    % Set coefficient based on direction of rotation
    Coe = strcmp(rightorleft, 'right') * 2 - 1;  % Coe = 1 for 'right', -1 for 'left'

    % Initialize the rotation parameters
    immproper = false;
    abc = [];

    % Check rotation matrix or axis-angle format and convert to Euler angles
    if isequal(size(rotm), [3 3])
        % Rotation matrix: 3x3
        if det(rotm) == -1
            % Inversion detected (improper rotation)
            immproper = true;
            rotm = -rotm; % Make the rotation proper
        end
        abc = Oper.Rotation2eul(rotm); % Convert to Euler angles (alpha, beta, gamma)
        
    elseif isequal(size(rotm), [4 1]) || isequal(size(rotm), [1 4])
        % Axis-angle format: 4x1 or 1x4
        if sym(rotm(end)) == -1
            immproper = true;
        end
        if sym(abs(rotm(end))) ~= 1
            % bug fix
            rotm(4) = -rotm(4);
            %
            abc = Oper.axang2eul(rotm(1:4)); % Convert axis-angle to Euler angles
            fprintf('Immproper input: Euler angle ZYZ format expected [alpha beta gamma]\n');
            immproper = true;
        else
            abc = rotm(1:3); % Axis vector in [nx, ny, nz] and angle in last element
        end
        
    elseif isequal(size(rotm), [5 1]) || isequal(size(rotm), [1 5])
        % Extended axis-angle (5x1 or 1x5)
        if sym(abs(rotm(end))) ~= 1
            error('Invalid input: Axis-angle format should be [nx ny nz theta det()]');
        end
        if sym(rotm(end)) == -1
            immproper = true;
        end
        % bug fix
        rotm(4) = -rotm(4);
        %
        abc = Oper.axang2eul(rotm(1:4)); % Convert axis-angle to Euler angles
    end

    % Define the inverse matrix based on whether the rotation is improper
    if immproper
        Invmat = ParityMat(SpinObj);  % Use parity matrix for improper rotation
    else
        Invmat = eye(size(SpinObj, 1));  % Identity matrix for proper rotation
    end

    % Use symbolic computation if requested
    if options.sym
        abc = sym(abc);
    end

    % Calculate the Wigner D matrix for the rotation
    U = WignerD(SpinObj, abc, rightorleft);

    % Apply the inverse matrix for improper rotations
    U = U * Invmat;

    % Simplify the result if symbolic computation is used
    if isa(U, 'sym')
        U = simplify(U);
    end
end
