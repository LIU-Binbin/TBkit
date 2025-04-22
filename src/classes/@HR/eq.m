function varargout = eq(H_hr1,H_hr2)
% EQ Compare two HR objects for equality
%
%   varargout = EQ(H_hr1,H_hr2) compares two HR objects and returns
%   logical result and potentially reshaped second Hamiltonian.
%
%   INPUT ARGUMENTS:
%       H_hr1, H_hr2 - HR objects to compare
%
%   OUTPUT ARGUMENTS (varargout):
%       1: logical_num - True if objects are equal
%       2: H_hr2 - Reshaped version of second Hamiltonian (if applicable)
%
%   NOTES:
%       - Compares NRPTS, WAN_NUM, vectorL, HcoeL and HnumL
%       - Returns reshaped H_hr2 when vectorL matches but is reordered
%
%   SEE ALSO:
%       HR, reseq
%

logical_num = true;
if H_hr1.NRPTS ~= H_hr2.NRPTS
    logical_num = false;
    varargout{1} = logical_num;
    return;
end
if H_hr1.WAN_NUM ~= H_hr2.WAN_NUM
    logical_num = false;
    varargout{1} = logical_num;
    return;
end
if H_hr1.NRPTS ~= H_hr2.NRPTS
    logical_num = false;
    varargout{1} = logical_num;
    return;
end
[issame_list,reshape_list] = ismember(H_hr1.vectorL,H_hr2.vectorL,'rows');
if sum(issame_list) ~= H_hr1.NRPTS
    logical_num = false;
    varargout{1} = logical_num;
    return;
else
    reshape_list = reshape_list';
    H_hr2 = H_hr2.reseq(':',reshape_list);
    varargout{2} = H_hr2 ;
end
if H_hr1.HcoeL ~= H_hr2.HcoeL
    logical_num = false;
    varargout{1} = logical_num;
    varargout{2} = (H_hr1.HcoeL == H_hr2.HcoeL);
    return;
end
if ~isequal(H_hr1.HnumL ,H_hr2.HnumL)
    logical_num = false;
    varargout{1} = logical_num;
    return;
end
varargout{1} = logical_num;
end
