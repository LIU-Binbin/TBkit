function C = directsum(A, B)
    % DIRECTSUM Direct sum (Kronecker sum) of two matrices.
    %   C = DIRECTSUM(A, B) computes the direct sum of matrices A and B.
    %
    %   Example:
    %       A = [1,2; 3,4];
    %       B = [5,6];
    %       C = directsum(A, B);
    
    [m, n] = size(A);
    [p, q] = size(B);
    
    % Create zero matrices with appropriate dimensions and class
    ZerosMat1 = zeros(m, q, class(A));
    ZerosMat2 = zeros(p, n, class(A));
    
    % Construct the direct sum matrix
    C = [A, ZerosMat1; 
          ZerosMat2, B];
end
% function C = directsum(A,B)
% SizeA = size(A);
% SizeB = size(B);
% ZerosMat = zeros(SizeA(1),SizeB(2),class(A));
% C = [A,ZerosMat ; ZerosMat.', B];
% end