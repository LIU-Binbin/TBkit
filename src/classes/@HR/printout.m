function Hout = printout(H_hr,print_list,mode)
% PRINTOUT Display or return HR object contents
%
%   HOUT = PRINTOUT(H_HR,PRINT_LIST,MODE) displays or returns HR contents
%
%   Inputs:
%       H_hr - HR object to display
%       print_list - Indices of terms to display [default: all]
%       mode - Display mode [default: 'HcoeL']
%   Output:
%       Hout - Requested matrix data
%
%   Notes:
%       - Modes: 'HnumL-sum', 'HnumL', 'HcoeL', 'HcoeL-sum'
%       - Can display individual terms or summed matrices
%       - Shows vector information with matrices
if nargin < 2
print_list = 1:H_hr.NRPTS;
end
if nargin < 3
mode = 'HcoeL';
end
switch mode
case 'HnumL-sum'
Hout = sum(H_hr.HnumL,3);
case 'HnumL'
for i = 1:length(print_list)
disp([H_hr.vectorL(print_list(i),:)]);
disp(H_hr.HnumL(:,:,print_list(i)));
end
Hout = H_hr.HnumL(:,:,print_list);
case 'HcoeL'
for i = 1:length(print_list)
disp([H_hr.vectorL(print_list(i),:)]);
disp(H_hr.HcoeL(:,:,print_list(i)));
end
Hout = H_hr.HcoeL(:,:,print_list);
case 'HcoeL-sum'
Hout =  sum(H_hr.HcoeL,3);
otherwise
fprintf("usage: H_hr.printout([1,3,4],'HcoeL')\n");
fprintf("usage: H_hr.printout([1,3,4],'HnumL')\n");
fprintf("Optional: HcoeL,HnumL,HnumL-sum,HcoeL-sum,\n");
end
end
