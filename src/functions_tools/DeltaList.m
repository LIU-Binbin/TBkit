function OutPutM = DeltaList(LeftL, RightL)
    % DeltaList constructs a binary matrix where the (i, j)-th entry is 1
    % if the i-th row of LeftL matches the j-th row of RightL, otherwise 0.
    %
    % Input:
    %   LeftL  - Matrix or cell array, with rows representing elements to match.
    %   RightL - Matrix or cell array, with rows representing elements to check against.
    %
    % Output:
    %   OutPutM - Binary matrix (Ncol x Ncol) where Ncol = number of rows in RightL.

    % Get the number of rows in RightL (Ncol)
    Ncol = size(RightL, 1);

    % Initialize output matrix with zeros
    OutPutM = zeros(Ncol);

    % Find rows in LeftL that match rows in RightL
    [ia, ib] = ismember(LeftL, RightL, 'rows');

    % ia gives logical indices of matches, ib gives corresponding indices in RightL
    iL = ia;  % indices of rows in LeftL that are matched
    jL = ib(ia);    % corresponding indices in RightL

    % Set the matching positions to 1 in the output matrix
    OutPutM(iL, jL) = 1;
end

% function OutPutM = DeltaList(LeftL,RightL)
%     % Width = size(LeftL,2);
%     Ncol = size(RightL,1);
%     OutPutM = zeros(Ncol);
%     [ia,ib] = ismember(LeftL,RightL,'rows');
%     iL = find(ia);
%     jL = ib(ia);
%     for i = 1:length(iL)
%         OutPutM(iL(i),jL(i)) = 1;
%     end
% end