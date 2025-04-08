function [outstring] = rm_line(filename, rm_number_list, outfilename)
    % RM_LINE: Removes specific lines from a file based on the provided line numbers.
    %   filename - The path to the input file.
    %   rm_number_list - A list of line numbers to be removed.
    %   outfilename - The output file where the modified content will be saved.
    %
    % The function will remove the specified lines from the input file and
    % save the result to outfilename. If outfilename is not provided, 
    % the result will not be written to a file.

    import linux_matlab.*;
    
    % Default mode is 'silence' (no output to console)
    if nargin < 3
        mode = 'silence';
    else
        mode = 'aloud';
    end
    
    % Read the file contents into a cell array of strings
    outstring = read_file(filename);
    
    % Remove the lines specified in rm_number_list
    outstring(rm_number_list) = [];
    
    % If mode is 'aloud', write the result to the output file
    if ~strcmp(mode, 'silence')
        write_file(outstring, outfilename);
    end
end
