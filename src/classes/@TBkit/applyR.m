function H_hk = applyR(H_hk, R)
%APPLYR Apply transformation matrix R to HK object coefficients
%   H_hk = APPLYR(H_hk, R) applies the transformation matrix R to both
%   numerical (HnumL) and symbolic (HcoeL) coefficient layers of an HK object.
%   The transformation is performed for each polynomial order in the HK object's
%   expansion using page-wise matrix multiplication.
%
% Input Arguments:
%   H_hk : HK object
%       Input Hamiltonian object containing coefficient matrices
%   R    : double | sym
%       Transformation matrix to apply (numeric or symbolic)
%
% Output Arguments:
%   H_hk : HK object
%       Transformed HK object with updated coefficients
%
% Methods Called:
%   HK.matgen          - Generates matrix cell array for transformation
%   HK.orderlist       - Gets polynomial order indices
%   HK.matrixtimespage - Performs page-wise matrix multiplication
%
% Example:
%   % Create HK object and identity transformation matrix
%   Hk = HK(...); 
%   R = eye(3);
%   % Apply transformation
%   transformed_Hk = applyR(Hk, R);
%
% See also:
%   HK, HK.matgen, HK.orderlist, HK.matrixtimespage

arguments
    H_hk HK        % HK object containing coefficient matrices
    R              % Transformation matrix (numeric/symbolic)
end

% Check numerical coefficients processing flag
if isequal(zeros(size(H_hk.HnumL)), H_hk.HnumL)
    num_label = false;  % No numerical coefficients
else
    num_label = true;   % Numerical coefficients exist
end

% Check symbolic coefficients processing flag
if isequal(sym(zeros(size(H_hk.HcoeL))), H_hk.HcoeL)
    coe_label = false;  % No symbolic coefficients
else
    coe_label = true;   % Symbolic coefficients exist
end

% Generate transformation matrix cell array
matcell = H_hk.matgen(R);

% Process numerical coefficients (HnumL)
if num_label
    for i = 0:H_hk.Degree
        Orderlist = HK.orderlist(i);  % Get order indices
        H_hk.HnumL(:,:,Orderlist) = HK.matrixtimespage(...
            matcell{i+1}, H_hk.HnumL(:,:,Orderlist));
    end
end

% Process symbolic coefficients (HcoeL)
if coe_label
    for i = 0:H_hk.Degree
        Orderlist = HK.orderlist(i);  % Get order indices
        H_hk.HcoeL(:,:,Orderlist) = HK.matrixtimespage(...
            matcell{i+1}, H_hk.HcoeL(:,:,Orderlist));
    end
end
end
