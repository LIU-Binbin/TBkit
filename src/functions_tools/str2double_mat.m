function [double_mat, str_mat] = str2double_mat(str)
    % str2double_mat - Converts a string matrix into a double matrix and a string matrix
    %
    % Syntax:
    %   [double_mat, str_mat] = str2double_mat(str)
    %
    % Input:
    %   str      - A string input that represents a matrix, with rows separated by newline characters ('\n') 
    %             and elements in each row separated by spaces or tabs.
    %
    % Output:
    %   double_mat - A matrix of doubles representing the input string matrix
    %   str_mat    - A matrix of strings corresponding to the input string matrix
    
    % Remove leading and trailing whitespace from the input string
    str = strtrim(str);

    % Split the input string into rows based on newline characters
    str_list = strsplit(str, '\n');
    
    % Initialize the double_mat matrix for conversion
    for i = 1:length(str_list)
        % Split each row string by spaces and convert to a double matrix
        double_mat(i, :) = str2double(strsplit(strtrim(str_list{i})));
    end
    
    % Convert the double matrix back to a string matrix
    str_mat = string(double_mat); % Convert to string type matrix
    
end
