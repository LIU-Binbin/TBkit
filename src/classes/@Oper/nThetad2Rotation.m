function RotationMat = nThetad2Rotation(thetad,n)
    if nargin < 2
        n = [0,0,1];
    end
    %             % Compute rotation matrices
    %             cth = cos(theta);
    %             sth = sin(theta);
    %             vth = (1 - cth);
    %             v = n;
    %             % Preallocate input vectors
    %             vx = zeros(1,1,numInputs,'like',axang);
    %             vy = vx;
    %             vz = vx;
    %
    %             % Shape input vectors in depth dimension
    %             vx(1,1,:) = v(:,1);
    %             vy(1,1,:) = v(:,2);
    %             vz(1,1,:) = v(:,3);
    %
    %             % Explicitly specify concatenation dimension
    %             tempR = cat(1, vx.*vx.*vth+cth,     vy.*vx.*vth-vz.*sth, vz.*vx.*vth+vy.*sth, ...
    %                 vx.*vy.*vth+vz.*sth, vy.*vy.*vth+cth,     vz.*vy.*vth-vx.*sth, ...
    %                 vx.*vz.*vth-vy.*sth, vy.*vz.*vth+vx.*sth, vz.*vz.*vth+cth);
    %
    %             R = reshape(tempR, [3, 3, length(vx)]);
    %             R = permute(R, [2 1 3]);
    n_x = n(1);
    n_y = n(2);
    n_z = n(3);
    %% left
    % RotationMat(1,1)=n_x^2*(1-cosd(theta))+cosd(theta);
    % RotationMat(1,2)=n_x*n_y*(1-cosd(theta))+n_z*sind(theta);
    % RotationMat(1,3)=n_x*n_z*(1-cosd(theta))-n_y*sind(theta);
    % RotationMat(2,1)=n_y*n_z*(1-cosd(theta))-n_z*sind(theta);
    % RotationMat(2,2)=n_y^2*(1-cosd(theta))+cosd(theta);
    % RotationMat(2,3)=n_y*n_z*(1-cosd(theta))+n_x*sind(theta);
    % RotationMat(3,1)=n_z*n_x*(1-cosd(theta))+n_y*sind(theta);
    % RotationMat(3,2)=n_z*n_y*(1-cosd(theta))-n_x*sind(theta);
    % RotationMat(3,3)=n_z^2*(1-cosd(theta))+cosd(theta);
    %% right
    RotationMat(1,1)=n_x^2*(1-cosd(thetad))+cosd(thetad);
    RotationMat(1,2)=n_x*n_y*(1-cosd(thetad))-n_z*sind(thetad);
    RotationMat(1,3)=n_x*n_z*(1-cosd(thetad))+n_y*sind(thetad);
    RotationMat(2,1)=n_y*n_x*(1-cosd(thetad))+n_z*sind(thetad);
    RotationMat(2,2)=n_y^2*(1-cosd(thetad))+cosd(thetad);
    RotationMat(2,3)=n_y*n_z*(1-cosd(thetad))-n_x*sind(thetad);
    RotationMat(3,1)=n_z*n_x*(1-cosd(thetad))-n_y*sind(thetad);
    RotationMat(3,2)=n_z*n_y*(1-cosd(thetad))+n_x*sind(thetad);
    RotationMat(3,3)=n_z^2*(1-cosd(thetad))+cosd(thetad);
end