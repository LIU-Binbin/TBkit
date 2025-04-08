function flag = strcontain(str, containlist)
    % Initialize the flag to 0 (false)
    flag = 0;

    % Iterate through each substring in the list
    for i = 1:length(containlist)
        % Check if the current substring is empty
        if isempty(containlist{i})
            % If the substring is empty, set flag to 1 (true) and exit the loop
            flag = 1;
            break;
        end
        % Check if the current substring is present in the main string
        if contains(str, containlist{i})
            % If found, set flag to 1 (true) and exit the loop
            flag = 1;
            break;
        end
    end
end
