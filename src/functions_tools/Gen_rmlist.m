function [orb_rm_index_list, orb_rm_list] = Gen_rmlist(orb_list, rm_or_not_functionname)
% Function to generate a list of orbitals to remove based on a user-defined judge function.
%
% Usage:
%   [orb_rm_index_list, orb_rm_list] = Gen_rmlist(orb_list, rm_or_not_functionname)
%   where `orb_list` is a list of orbitals, and `rm_or_not_functionname` is the function
%   used to determine which orbitals to remove.
%
% The judge function should take a single orbital `rc` as input and return 1 (remove) or 0 (do not remove).
%
% Inputs:
%   orb_list               - List of orbitals (each row is an orbital with its coordinates).
%   rm_or_not_functionname - Name of the user-defined function to judge if an orbital should be removed (default: 'rm_or_not').
%
% Outputs:
%   orb_rm_index_list - Indices of orbitals to remove.
%   orb_rm_list       - List of orbitals to remove.

    %% Check the number of input arguments
    if nargin < 2
        rm_or_not_functionname = 'rm_or_not';  % Default judge function
        % Optionally, provide a template for the judge function here.
    end
    
    %% Check if the specified judge function exists
    rm_or_not_functionname_m = rm_or_not_functionname + ".m";
    if ~exist(rm_or_not_functionname_m, 'file')
        error('The file for the specified judge function does not exist! Please provide the correct rm_or_not function.');
    end
    
    %% Test if the function works correctly
    disp('Testing if the function works with the default input rc = [0, 0, 0]...');
    rc_test = [0, 0, 0];
    test = rm_or_not_eval(rm_or_not_functionname, rc_test);
    
    % Check if the function returns 0 or 1
    if test == 0 || test == 1
        disp('Test passed!');
    else
        disp('Test failed.');
        error('The rm_or_not function does not work properly. Please check and try again.');
    end
    
    %% Generate orb_rm_index_list and orb_rm_list
    orb_rm_index_list = [];  % Initialize index list
    orb_rm_list = [];        % Initialize the list of orbitals to remove
    
    [norb, ~] = size(orb_list);  % Get the number of orbitals
    for i = 1:norb
        temp_orb = orb_list(i, :);  % Extract each orbital
        logic_result = rm_or_not_eval(rm_or_not_functionname, temp_orb);  % Apply the judge function
        
        % If the orbital should be removed (i.e., logic_result == 1)
        if logic_result == 1
            orb_rm_index_list = [orb_rm_index_list; i];  % Add index to removal list
            orb_rm_list = [orb_rm_list; temp_orb];        % Add orbital to removal list
        end
    end
    
    %% Return the lists
end

%% Helper function to evaluate the judge function
% This function dynamically evaluates the user-defined judge function with the given orbital (rc).
function logic_result = rm_or_not_eval(rm_or_not_functionname, rc)
    % Construct the command string to call the user-defined function with the input orbital (rc)
    eval_string = rm_or_not_functionname + "(" + num2str(rc(1)) + "," + num2str(rc(2)) + "," + num2str(rc(3)) + ");";
    
    % Evaluate the function and get the result (0 or 1)
    logic_result = eval(eval_string);
end
