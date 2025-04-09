function StrM = str2mat_overwrited(StrL, delimiter, dim)
% STR2MAT Converts a string array into a matrix by splitting each string.
% Each string is split by the specified delimiter, and the result is stored
% in a matrix where each row corresponds to a string and each column contains
% a part of the split string.
%
% Input:
%   StrL      - A string array where each row is a string to be split.
%   delimiter - The delimiter used to split the strings (default: ' ').
%   dim       - Dimension to split along (default: 2, meaning split by columns).
%
% Output:
%   StrM      - A matrix where each row contains the split parts of the input strings.

    % Set default values for optional arguments
    if nargin < 2
        delimiter = ' ';
    end
    if nargin < 3
        dim = 2;
    end

    % Number of input strings
    nStrL = size(StrL, 1);

    % Initialize variables to store split results and maximum length
    tmpStrL = cell(nStrL, 1); % Temporary storage for split strings
    checklength = zeros(nStrL, 1); % Store the length of each split result
    maxlength = 0; % The maximum length of split parts

    % Split each string and calculate the maximum length
    for i = 1:nStrL
        tmpStrL{i} = split(strtrim(StrL(i, :)), delimiter, dim);
        checklength(i) = length(tmpStrL{i});
        maxlength = max(maxlength, checklength(i));
    end

    % Create an empty matrix with the maximum length
    StrM = repmat("", [nStrL, maxlength]);

    % Fill the matrix with the split parts
    for i = 1:nStrL
        StrM(i, 1:checklength(i)) = tmpStrL{i};
    end

    % Call external function to clean the matrix
    StrM = cleanStrM(StrM);
end



% 示例输入
% StrL = ["Hello World MATLAB"; "This is a test"; "Another test case"];
% 
% % 调用函数，默认使用空格分隔
% StrM = str2mat(StrL, ' ', 2);
% 
% disp(StrM);



% function StrM = str2mat(StrL,delimiter,dim)
% arguments
%     StrL string;
%     delimiter  = ' ';
%     dim =2;
% end
% maxlength = 1;
% nStrL = size(StrL,1);
% for i = 1:nStrL
%     %
%     tmpStrL{i} = split(strtrim(StrL(i,:)),delimiter,dim);
%     checklength{i} = length(tmpStrL{i});
%     if checklength{i} > maxlength
%         maxlength = checklength{i};
%     else
% 
%     end
% end
% StrM = repmat("",[nStrL,maxlength]);
% for i = 1:nStrL
%     StrM(i,1:checklength{i}) = tmpStrL{i};
% end
% StrM = park.cleanStrM(StrM);
% end