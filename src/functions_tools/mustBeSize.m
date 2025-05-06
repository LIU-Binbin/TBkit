function mustBeSize(a, b)
    % MUSTBESIZE validates the size of input a against the expected size b.
    % This function checks if the size of a matches the given size b.
    
    % Case 1: If b is a 1xN vector, check if size(a) matches size(b)
    if size(b, 1) == 1
        if ~isequal(size(a), b)
            eid = 'Size:notRequired';
            msg = ['Inputs must have size ' mat2str(b) ', but size(a) is ' mat2str(size(a))];
            throwAsCaller(MException(eid, msg));
        end
        
    % Case 2: If b is a Nx2 matrix, check if size(a) matches either of the rows of b
    elseif size(b, 2) == 2
        Check = false;
        for i = 1:size(b,1)
            if isequal(size(a), b(i,:))
                Check = true;
            end
        end
        if ~Check
            eid = 'Size:notRequired';
            msg = ['Inputs must have size ' mat2str( b(i,:)) ', but size(a) is ' mat2str(size(a))];
            throwAsCaller(MException(eid, msg));
        end
    else
        % Case when b is neither 1xN nor 2xN, handle accordingly
        % (you can add more conditions if needed)
        eid = 'Size:notSupported';
        msg = 'Input b must be either a 1xN or 2xN matrix.';
        throwAsCaller(MException(eid, msg));
    end
end
