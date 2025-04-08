function TEST_example_in_TBkit(dir, runfile, basedir)
    % TEST_example_in_TBkit - A function to set up and run an example in the VASP library
    %
    % Syntax:
    %   TEST_example_in_TBkit(dir, runfile, basedir)
    %
    % Input:
    %   dir      - Directory name (relative to 'basedir') containing the example files
    %   runfile  - The script or command to run after setting up the example
    %   basedir  - The base directory from where to copy example files (default is '../example/')
    %
    % Description:
    %   This function copies the specified example directory into a test folder, 
    %   runs a specified script, and cleans up the test folder afterwards.
    
    % Store current working directory to return after the test
    workshop = pwd;
    
    % Create a 'test' directory for the test run
    mkdir('test');
    cd test;
    
    % If basedir is not provided, default to '../example/'
    if nargin < 3
        basedir = '../example/';
    end
    
    % Copy the example directory to the test directory
    copyfile([basedir, dir]);
    
    % If runfile is not provided, display a message to check the runfile
    if nargin < 2
        fprintf('No runfile provided. Please specify a valid runfile.\n');
        return;
    end
    
    % Try running the specified runfile
    try
        run(runfile);
    catch
        % Display an error message if the runfile fails to execute
        fprintf('An error occurred while running the script. Please contact me: parkman@buaa.edu.cn\n');
    end
    
    % Return to the original directory
    cd(workshop);
    
    % Clean up by removing the 'test' directory
    rmdir('test', 's');
end
