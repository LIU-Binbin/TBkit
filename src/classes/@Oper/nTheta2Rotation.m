function RotationMat = nTheta2Rotation(theta, n)
    % nTheta2Rotation Convert axis-angle rotation representation to rotation matrix
    %   This function converts a 3D rotation specified in axis-angle form (theta, n)
    %   to an orthonormal rotation matrix. The input theta is a scalar representing
    %   the rotation angle in radians, and n is a 3-element vector specifying
    %   the rotation axis. The output is a 3-by-3 orthonormal rotation matrix.
    %
    %   Syntax:
    %       RotationMat = nTheta2Rotation(theta, n)
    %       RotationMat = nTheta2Rotation(theta)  % Uses default n = [0,0,1]
    %
    %   Inputs:
    %       theta - Rotation angle in radians (scalar)
    %       n     - Rotation axis (3-element vector) [default = [0,0,1]]
    %
    %   Outputs:
    %       RotationMat - 3-by-3 orthonormal rotation matrix
    %
    %   Example:
    %       % Convert a rotation from axis-angle to rotation matrix
    %       axang = [0 1 0 pi/2];
    %       R = nTheta2Rotation(axang(4), axang(1:3));
    
    % Set default rotation axis if not provided
    if nargin < 2
        n = [0, 0, 1];  % Default rotation axis: z-axis
    end
    
    % Normalize the rotation axis
    n = n./norm(n);  % Ensure the axis vector is unit-length
    
    % Convert rotation angle from radians to degrees
    thetaDegrees = rad2deg(theta);  % Convert theta from radians to degrees
    
    % Compute the rotation matrix using the axis-angle representation
    RotationMat = Oper.nThetad2Rotation(thetaDegrees, n);  
end