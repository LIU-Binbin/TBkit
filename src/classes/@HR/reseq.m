function H_hr = reseq(H_hr,wan_list,nrpt_list,options)
% RESEQ Reorder Wannier functions and hopping terms in HR object
%
%   H_HR = RESEQ(H_HR,WAN_LIST,NRPT_LIST) reorders Wannier functions
%   and hopping terms in the HR object
%
%   Inputs:
%       H_hr - HR object to reorder
%       wan_list - List of Wannier indices to reorder [default: ':']
%       nrpt_list - List of hopping terms to reorder [default: ':']
%   Output:
%       H_hr - Reordered HR object
%
%   Notes:
%       - Maintains consistency across all object properties
%       - Handles 'sparse', 'list', and 'mat' storage types
%       - Preserves overlap terms when present
arguments
    H_hr  HR;
    wan_list  = ':';
    nrpt_list = ':';
    options.OverlapTest = true;
end
    
if H_hr(1).overlap && options.OverlapTest
    H_hr(1) = add_empty_one(H_hr(1),wan_list,nrpt_list,'OverlapTest',false);
    H_hr(2) = add_empty_one(H_hr(2),wan_list,nrpt_list,'OverlapTest',false);
    return
end



if ~isequal(wan_list,':')
    if strcmp(H_hr.Type,'sparse')
        for i = 1:H_hr.NRPTS
            H_hr.HnumL{i}=H_hr.HnumL{i}(wan_list,wan_list);
        end
    elseif strcmp(H_hr.Type,'list')
        wan_list =  all(ismember(H_hr.vectorL(:,H_hr.Dim+1:H_hr.Dim+2),(wan_list)),2);
        if H_hr.num
            H_hr.HnumL=H_hr.HnumL(wan_list,:);
        end
        if H_hr.coe
            H_hr.HcoeL=H_hr.HcoeL(wan_list,:);
        end
        H_hr.vectorL=H_hr.vectorL(wan_list,:);
        if ~isempty(H_hr.vectorL_map)
            H_hr.vectorL_map=H_hr.vectorL_map(wan_list,:);
        end
    elseif strcmp(H_hr.Type,'mat')
        if H_hr.num
            H_hr.HnumL=H_hr.HnumL(wan_list,wan_list,:);
        else
            H_hr.HnumL = [];
            if H_hr.overlap
                H_hr.SnumL=[];
            end
        end
        if H_hr.coe
            H_hr.HcoeL=H_hr.HcoeL(wan_list,wan_list,:);
        else
            H_hr.HcoeL = sym([]);
        end
    end
    if ~isempty(H_hr.sites)
        try
            H_hr.sites = H_hr.sites(wan_list);
        catch
        end
    end
    if ~isempty( H_hr.orbL )
        H_hr.orbL = H_hr.orbL(wan_list,:);
    end
    if ~isempty( H_hr.elementL )
        H_hr.elementL = H_hr.elementL(wan_list,:);
    end
    if ~isempty( H_hr.quantumL)
        H_hr.quantumL = H_hr.quantumL(wan_list,:);
    end
end
if ~isequal(nrpt_list,':')
    H_hr.vectorL=H_hr.vectorL(nrpt_list,:);
    if ~isempty(H_hr.vectorL_map)
        try
            rows = num2cell(H_hr.vectorL, 2);
            % 对每一行应用 mat2str 函数
            H_hr.vectorL_map  = cellfun(@mat2str, rows, 'UniformOutput', false);
            H_hr.vectorL_map=H_hr.vectorL_map(nrpt_list,:);
        catch
        end
    end
    if strcmp(H_hr.Type,'sparse')
        H_hr.HnumL=H_hr.HnumL(nrpt_list);
    elseif strcmp(H_hr.Type,'list')
        if H_hr.vectorhopping
            H_hr.AvectorL=H_hr.AvectorL(nrpt_list,:);
            H_hr.BvectorL=H_hr.BvectorL(nrpt_list,:);
            CL1 = H_hr.CvectorL(1:end/2,:);
            CL2 = H_hr.CvectorL(end/2+1:end,:);
            H_hr.CvectorL = [CL1(nrpt_list,:);CL2(nrpt_list,:)];
            return;
        end
        if H_hr.num
            H_hr.HnumL=H_hr.HnumL(nrpt_list);
        end
        if H_hr.coe
            H_hr.HcoeL=H_hr.HcoeL(nrpt_list);
        end
    elseif strcmp(H_hr.Type,'mat')
        if H_hr.num
            H_hr.HnumL=H_hr.HnumL(:,:,nrpt_list);
        end
        if H_hr.coe
            H_hr.HcoeL=H_hr.HcoeL(:,:,nrpt_list);
        end
    end
end
end
