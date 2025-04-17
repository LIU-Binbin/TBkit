function H_hr = set_hop(H_hr, amp, hi, hj, vector_list, mode)
%SET_HOP Updates hopping terms in Hamiltonian data structure
%   H_hr = SET_HOP(H_hr, amp, hi, hj, vector_list, mode) modifies hopping terms:
%       amp: Hopping amplitude(s) - scalar, matrix or spinor components
%       hi: Initial Wannier function indices (row/column indices)
%       hj: Final Wannier function indices (row/column indices)
%       vector_list: Lattice vectors [nx, ny, nz] for hopping terms
%       mode: Operation mode ('set', 'add', 'sym', 'symadd')
%
%   Supports different storage types: matrix, sparse, list
%   Handles spinor components (2x2, 4x4, 8x8 matrices) via amplitude dimensions

    %% Input validation and initialization
    arguments
        H_hr;          % Hamiltonian data structure
        amp;           % Hopping amplitude(s)
        hi;            % Initial state indices
        hj;            % Final state indices
        vector_list;   % Lattice vectors [nx, ny, nz]
        mode = 'set';  % Operation mode (default: overwrite)
    end
    
    numWannHalf = H_hr.WAN_NUM/2;         % Half the Wannier count for spin systems
    [nVectors, ~] = size(vector_list);    % Number of lattice vectors
    
    %% Configure symmetry handling
    if contains(mode, 'sym')
        H_hr.num = false;   % Disable numerical storage
        H_hr.coe = true;    % Enable symbolic coefficients
    end
    
    %% Process hopping terms for each lattice vector
    for iVec = 1:nVectors
        currentVector = vector_list(iVec, :);
        
        %% Matrix-based processing
        if isMatrixOperation(hi, hj, H_hr)
            H_hr = processMatrixForm(currentVector, H_hr, mode, amp, hi, hj);
            
        %% Sparse matrix handling    
        elseif isSparseOperation(hi, hj, H_hr)
            H_hr = processSparseForm(currentVector, H_hr, mode, amp, hi, hj);
            
        %% List-based processing
        elseif isListOperation(hi, hj, H_hr)
            H_hr = processListForm(currentVector, H_hr, mode, amp, hi, hj);
            
        %% Specialized spin/band handling
        else
            H_hr = handleSpecialCases(currentVector, H_hr, mode, amp, hi, hj, numWannHalf);
        end
    end
    
    %% Nested helper functions %%
    function flag = isMatrixOperation(hi, hj, H_hr)
        flag = length(hi) == length(hj) && length(hi) > 1 && strcmp(H_hr.Type, 'mat');
    end
    
    function H_hr = processMatrixForm(vector, H_hr, mode, amp, hi, hj)
        switch mode
            case {'set', 'add'}
                ampMat = zeros(H_hr.WAN_NUM);
            case {'sym', 'symadd'}
                ampMat = sym(zeros(H_hr.WAN_NUM));
        end
        ampMat(sub2ind(size(ampMat), hi, hj)) = amp;
        H_hr = H_hr.set_hop_mat(ampMat, vector, mode);
    end
    
    function flag = isSparseOperation(hi, hj, H_hr)
        flag = length(hi) == length(hj) && strcmp(H_hr.Type, 'sparse');
    end
    
    function H_hr = processSparseForm(vector, H_hr, mode, amp, hi, hj)
        ampMat = sparse(hi, hj, amp, H_hr.WAN_NUM, H_hr.WAN_NUM);
        H_hr = H_hr.set_hop_mat(ampMat, vector, mode);
    end
    
    function flag = isListOperation(hi, hj, H_hr)
        flag = length(hi) == length(hj) && length(hi) > 1 && strcmp(H_hr.Type, 'list');
    end
    
    function H_hr = processListForm(vector, H_hr, mode, amp, hi, hj)
        for idx = 1:length(hi)
            H_hr = H_hr.set_hop_single(amp(idx), hi(idx), hj(idx), vector, mode);
        end
    end
    
    function H_hr = handleSpecialCases(vector, H_hr, mode, amp, hi, hj, numWannHalf)
        switch size(amp, 1)
            case 1
                %                     % single mode
                H_hr = H_hr.set_hop_single(amp,hi,hj,vector,mode);
            case 2  % Spin-1/2 system
                H_hr = updateSpinComponents(vector, H_hr, mode, amp, hi, hj, numWannHalf);
            case 4  % Four-band system
                warning('Four-band mode not implemented');
            case 8  % Eight-band system
                warning('Eight-band mode not implemented');
            otherwise  % Generic matrix input
                H_hr = H_hr.set_hop_mat(amp, vector, mode);
        end
    end
    
    function H_hr = updateSpinComponents(vector, H_hr, mode, amp, hi, hj, numWannHalf)
        % Update spin-up components
        H_hr = H_hr.set_hop_single(amp(1,1), hi, hj, vector, mode);
        H_hr = H_hr.set_hop_single(amp(1,2), hi, hj+numWannHalf, vector, mode);
        % Update spin-down components
        H_hr = H_hr.set_hop_single(amp(2,1), hi+numWannHalf, hj, vector, mode);
        H_hr = H_hr.set_hop_single(amp(2,2), hi+numWannHalf, hj+numWannHalf, vector, mode);
    end
end