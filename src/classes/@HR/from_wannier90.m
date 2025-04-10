function H_hr = from_wannier90(filename,Type,options)
% FROM_WANNIER90 Create HR object from Wannier90 output
%
%   H_hr = FROM_WANNIER90(filename,Type,options) imports Wannier90
%   _hr.dat file into an HR object.
%
%   INPUT ARGUMENTS:
%       filename - Path to wannier90_hr.dat file
%       Type - Output type ('mat' or 'list')
%       options - Structure with parameters:
%           Accuracy: Numerical threshold for filtering
%           overlap: Include overlap matrix (logical)
%
%   OUTPUT ARGUMENTS:
%       H_hr - HR object constructed from Wannier90 data
%
%   NOTES:
%       - Handles both regular and overlap Hamiltonians
%       - Can apply degeneracy factors from NRPT_list
%       - Filters small elements based on Accuracy
%
%   SEE ALSO:
%       HR, hrdat_read
%
%   AUTHOR:
%       [Your Name] ([Your Email])
%       [Creation Date]

arguments
filename  = 'wannier90_hr.dat';
Type  char = 'mat';
options.Accuracy = 1e-6;
options.overlap = false
end
if options.overlap
[dataArray,NRPT_list,NRPTS,NUM_WAN]=HR.hrdat_read(filename{1});
[dataArray2,NRPT_list_S,~,~]=HR.hrdat_read(filename{2});
if strcmp(Type ,'mat')
Vec_Fir = dataArray{:, 1};
Vec_Sec = dataArray{:, 2};
Vec_Thi = dataArray{:, 3};
h_real = dataArray{:, 6};
h_imag = dataArray{:, 7};
V_f = reshape(Vec_Fir, NUM_WAN*NUM_WAN, NRPTS);
V_s = reshape(Vec_Sec, NUM_WAN*NUM_WAN, NRPTS);
V_t = reshape(Vec_Thi, NUM_WAN*NUM_WAN, NRPTS);
vectorL=[V_f(1,:)',V_s(1,:)',V_t(1,:)'];
HnumL_real = reshape(h_real, NUM_WAN,NUM_WAN, NRPTS);
HnumL_imag = reshape(h_imag, NUM_WAN,NUM_WAN, NRPTS);
HnumL = HnumL_real+1i*HnumL_imag;
if strcmp(Type,'n')
for i = 1 : NRPTS
HnumL(:,:,i) = HnumL(:,:,i)/NRPT_list(i);
end
else
end
Vec_Fir = dataArray2{:, 1};
Vec_Sec = dataArray2{:, 2};
Vec_Thi = dataArray2{:, 3};
h_real = dataArray2{:, 6};
h_imag = dataArray2{:, 7};
V_f = reshape(Vec_Fir, NUM_WAN*NUM_WAN, NRPTS);
V_s = reshape(Vec_Sec, NUM_WAN*NUM_WAN, NRPTS);
V_t = reshape(Vec_Thi, NUM_WAN*NUM_WAN, NRPTS);
vectorL_overlap=[V_f(1,:)',V_s(1,:)',V_t(1,:)'];
SnumL_real = reshape(h_real, NUM_WAN,NUM_WAN, NRPTS);
SnumL_imag = reshape(h_imag, NUM_WAN,NUM_WAN, NRPTS);
SnumL = SnumL_real+1i*SnumL_imag;
if strcmp(Type,'n')
for i = 1 : NRPTS
SnumL(:,:,i) = SnumL(:,:,i)/NRPT_list(i);
end
else
end
HcoeL = sym([]);
ScoeL = sym([]);
elseif strcmp(Type ,'list')
DATAARRAY = cell2mat(dataArray(1:7));
HOPARRAY = DATAARRAY(:,6)+1i*DATAARRAY(:,7);
HnumL_select = abs(HOPARRAY)>options.Accuracy;
vectorL = DATAARRAY(HnumL_select,1:H_hr.Dim+2);
[~,~,original_label] = unique(vectorL(:,1:H_hr.Dim),'rows');
HnumL = HOPARRAY(HnumL_select)./NRPT_list(original_label);
HcoeL = sym([]);
DATAARRAY = cell2mat(dataArray2(1:7));
OVERLAPARRAY = DATAARRAY(:,6)+1i*DATAARRAY(:,7);
SnumL_select = abs(OVERLAPARRAY)>options.Accuracy;
vectorL_overlap = DATAARRAY(SnumL_select,1:H_hr.Dim+2);
[~,~,original_label] = unique(vectorL_overlap(:,1:H_hr.Dim),'rows');
SnumL = OVERLAPARRAY(SnumL_select)./NRPT_list_S(original_label);
ScoeL = sym([]);
end
H_hr = HR(NUM_WAN,vectorL,'HnumL',HnumL,'HcoeL',HcoeL,'SnumL',SnumL,'ScoeL',ScoeL,'Type',Type,'overlap',true,'sym','false');
H_hr.num = true;
H_hr.coe = false;
H_hr.vectorL_overlap = vectorL_overlap;
H_hr.overlap = true;
else
[dataArray,NRPT_list,NRPTS,NUM_WAN]=HR.hrdat_read(filename);
if strcmp(Type ,'mat')
Vec_Fir = dataArray{:, 1};
Vec_Sec = dataArray{:, 2};
Vec_Thi = dataArray{:, 3};
h_real = dataArray{:, 6};
h_imag = dataArray{:, 7};
V_f = reshape(Vec_Fir, NUM_WAN*NUM_WAN, NRPTS);
V_s = reshape(Vec_Sec, NUM_WAN*NUM_WAN, NRPTS);
V_t = reshape(Vec_Thi, NUM_WAN*NUM_WAN, NRPTS);
vectorL=[V_f(1,:)',V_s(1,:)',V_t(1,:)'];
HnumL_real = reshape(h_real, NUM_WAN,NUM_WAN, NRPTS);
HnumL_imag = reshape(h_imag, NUM_WAN,NUM_WAN, NRPTS);
HnumL = HnumL_real+1i*HnumL_imag;
if strcmp(Type,'n')
for i = 1 : NRPTS
HnumL(:,:,i) = HnumL(:,:,i)/NRPT_list(i);
end
else
end
HcoeL = sym([]);
elseif strcmp(Type ,'list')
DATAARRAY = cell2mat(dataArray(1:7));
HOPARRAY = DATAARRAY(:,6)+1i*DATAARRAY(:,7);
HnumL_select = abs(HOPARRAY)>options.Accuracy;
vectorL = DATAARRAY(HnumL_select,1:H_hr.Dim+2);
[~,~,original_label] = unique(vectorL(:,1:3),'rows');
HnumL = HOPARRAY(HnumL_select)./NRPT_list(original_label);
HcoeL = sym([]);
end
H_hr = HR(NUM_WAN,vectorL,'HnumL',HnumL,'HcoeL',HcoeL,'Type',Type,'sym',false);
H_hr.num = true;H_hr.coe = false;
end
H_hr.Basis_num = H_hr.WAN_NUM;
end
