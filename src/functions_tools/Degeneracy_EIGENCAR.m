function [DEGENCAR, NODEINSECT] = Degeneracy_EIGENCAR(EIGENCAR, highK, dE)
% DEGGENERACY_EIGENCAR calculates the degeneracy and node intersections of eigenvalues.
% 
% Input:
%   EIGENCAR  - Matrix of eigenvalues (NBANDS x NhighK).
%   highK     - Indices of high-symmetry k-points.
%   dE        - Energy difference threshold for degeneracy checking (default: 0.1).
% 
% Output:
%   DEGENCAR  - Matrix of degeneracies (NBANDS x NhighK).
%   NODEINSECT - Matrix of node intersections between adjacent bands (NBANDS-1 x NhighK-1).

    % Default value for dE if not provided
    if nargin < 3
        dE = 0.1;
    end

    % Get the size of the eigenvalue matrix
    [NBANDS, ~] = size(EIGENCAR);
    NhighK = length(highK);

    % Initialize output matrices
    DEGENCAR = zeros(NBANDS, NhighK); 
    NODEINSECT = zeros(NBANDS-1, NhighK-1);

    % Calculate degeneracies for each k-point
    for i = 1:NhighK
        DEGENCAR(1, i) = 1; % The first band is always degenerate by itself
        count = 1;
        for j = 2:NBANDS
            if (EIGENCAR(j, highK(i)) - EIGENCAR(j-1, highK(i))) < dE
                DEGENCAR(count, i) = DEGENCAR(count, i) + 1; % Increase degeneracy count
            else
                count = count + 1; % Start counting a new degeneracy group
                DEGENCAR(count, i) = 1;
            end 
        end
    end

    % Calculate node intersections between adjacent bands at high symmetry k-points
    for i = 1:NBANDS-1
        for j = 1:NhighK-1
            % Check the absolute difference between adjacent eigenvalues at highK points
            NODEINSECT(i, j) = sum(abs(EIGENCAR(i+1, highK(j)+1:highK(j+1)-1) - ...
                                       EIGENCAR(i, highK(j):highK(j+1)-1)) < dE);
        end
    end
end


% function [DEGENCAR, NODEINSECT]= Degeneracy_EIGENCAR(EIGENCAR,highK,dE)
% if nargin <3
%     dE =0.1;
% end
%     [NBANDS,~] = size(EIGENCAR);
%     NhighK = length(highK);
%     DEGENCAR = zeros(NBANDS,NhighK);
%     NODEINSECT = zeros(NBANDS-1,NhighK-1);
%     % NhighK
%     for i = highK
%         DEGENCAR(1,i) = 1;
%         count = 1;
%         for j = 2:NBANDS
%             if (EIGENCAR(j,i)-EIGENCAR(j-1,i)) < dE
%                 DEGENCAR(count,i) = DEGENCAR(count,i)+1;
%             else
%                 count = count+1;
%                 DEGENCAR(count,i) = 1;
%             end 
%         end
%     end
%     % 
%     for i = 1:NBANDS-1
%         for j = 1:NhighK-1
%             NODEINSECT(i,j) = sum(abs(EIGENCAR(i+1,highK(j)+1:highK(j+1)-1)-...
%                 EIGENCAR(i,highK(j):highK(j+1)-1)) < dE);
%         end
%     end
% end