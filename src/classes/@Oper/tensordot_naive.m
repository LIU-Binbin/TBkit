function M = tensordot_naive(A,B,sizeLast)
if nargin < 3
sizeA = size(A);
sizeB = size(B);
sizeLast = sizeB(end);
end
if length(sizeA) ~= length(sizeB)
A = reshape(A,sqrt(length(A(:))/sizeLast),sqrt(length(A(:))/sizeLast),sizeLast);
end
M = A(:,:,1)* B(:,:,1);
for i = 2:sizeLast
M = M + A(:,:,i)* B(:,:,i);
end
end
