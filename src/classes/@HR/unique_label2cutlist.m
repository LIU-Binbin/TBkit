function cutlist= unique_label2cutlist(unique_label,NRPTS)
%UNIQUE_LABEL2CUTLIST Convert unique labels to cut list ranges
%
%   CUTLIST = UNIQUE_LABEL2CUTLIST(UNIQUE_LABEL, NRPTS)
%   Creates index ranges for groups of identical labels
%
%   Inputs:
%       unique_label - Array of unique group labels
%       NRPTS        - Total number of points
%
%   Output:
%       cutlist - N×2 array of [start, end] indices for each group
%
%   Note:
%       Last group always ends at NRPTS
cutlist(:,1)= unique_label;
cutlist(1:end-1,2)= unique_label(2:end)-1;
cutlist(end,2) = NRPTS;
end
