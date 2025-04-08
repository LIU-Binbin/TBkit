function H_hr = autohermi(H_hr,mode,options)
arguments
H_hr HR;
mode =  'sym';
options.enforce_list = true;
end
if nargin<2
mode = 'sym';
end
if options.enforce_list
if ~strcmp(H_hr.type, 'list')
H_hr = H_hr.rewrite();
giveback  = true;
else
giveback  = false;
end
else
giveback  = true;
end
H_hr_tmp = H_hr;
switch mode
case 'sym'
fprintf('For sym Hr, The requirement is strick!');
fprintf('The sym vasiable is real or ');
fprintf('in the original Hr, use conj()\n');
switch H_hr.type
case 'mat'
for i = 1:H_hr.NRPTS
vector_tmp = H_hr.vectorL(i,:);
vector_tmp_oppo = -vector_tmp;
[~,j]=ismember(vector_tmp_oppo,H_hr_tmp.vectorL,'rows');
fprintf('check %3d th : NRPT,\n the vector is %3d %3d %3d,\n the opposite vector is %3d,%3d,%3d the %3d th NRPT\n',i,...
vector_tmp(1),vector_tmp(2),vector_tmp(3),...
vector_tmp_oppo(1),vector_tmp_oppo(2),vector_tmp_oppo(3),j);
if i == j
if ~isequal(H_hr.HcoeL(:,:,i) ,H_hr_tmp.HcoeL(:,:,j)')
fprintf('The homecell hamilton is not hermi, never mind, we will hermi it enforcely!\n');
disp(H_hr.HcoeL(:,:,i));
H_hr_tmp = H_hr_tmp.set_hop_mat(H_hr.HcoeL(:,:,i)'/2+H_hr_tmp.HcoeL(:,:,j)/2,[0,0,0],'sym');
fprintf('change into\n');
disp(H_hr_tmp.HcoeL(:,:,i));
end
continue;
end
if j == 0
fprintf('The opposite vector hamilton does not exist, build it!\n');
H_hr_tmp = H_hr_tmp.set_hop_mat(H_hr.HcoeL(:,:,i)',-vector_tmp,'sym');
disp(H_hr_tmp.HcoeL(:,:,i)');
continue;
elseif ~isequal(H_hr.HcoeL(:,:,i) ,H_hr_tmp.HcoeL(:,:,j)')
fprintf('The opposite vector hamilton is not hermi, replace it by strong mat! \n');
N1 = nnz(H_hr.HcoeL(:,:,i));
N2 = nnz(H_hr_tmp.HcoeL(:,:,j)');
if N1 >= N2
fprintf('The %3d th NRPT is stonger!\n',i);
disp(H_hr.HcoeL(:,:,i));
H_hr_tmp = H_hr_tmp.set_hop_mat(H_hr.HcoeL(:,:,i)',-vector_tmp,'sym');
elseif N1 < N2
fprintf('The %3d th NRPT is stonger!\n',j);
disp(H_hr_tmp.HcoeL(:,:,j)');
H_hr_tmp = H_hr_tmp.set_hop_mat(H_hr_tmp.HcoeL(:,:,j)',vector_tmp,'sym');
else
end
elseif isequal(H_hr.HcoeL(:,:,i),H_hr_tmp.HcoeL(:,:,j)')
disp('hermi test pasts');
else
warning('!!!!!');
end
end
case 'list'
DIM = H_hr.Dim;
for i = 1:H_hr.NRPTS
vector_tmp = H_hr.vectorL(i,:);
vector_tmp_oppo(1:DIM) = -vector_tmp(1:DIM);
vector_tmp_oppo(DIM+1) = vector_tmp(DIM+2);
vector_tmp_oppo(DIM+2) = vector_tmp(DIM+1);
[~,j]=ismember(vector_tmp_oppo,H_hr_tmp.vectorL,'rows');
if i == j
if ~isequal(H_hr.HcoeL(i) ,H_hr_tmp.HcoeL(j)')
if H_hr.HcoeL(i)== sym(0)
FORMAT = strcat("The testing homecell hamilton ",fold(@strcat,repmat("
Input = num2cell(vector_tmp);
Input{DIM+3} = string(H_hr.HcoeL(j)');
fprintf(FORMAT,Input{:});
H_hr_tmp = H_hr_tmp.set_hop(H_hr.HcoeL(j)',vector_tmp(DIM+1),vector_tmp(DIM+2),vector_tmp(1:H_hr.Dim),'sym');
elseif H_hr.HcoeL(j)== sym(0)
FORMAT = strcat("The opposite homecell  hamilton ",fold(@strcat,repmat("
Input = num2cell(vector_tmp);
Input{DIM+3} = string(H_hr.HcoeL(i)');
fprintf(FORMAT, Input{:});
H_hr_tmp = H_hr_tmp.set_hop(H_hr.HcoeL(i)',vector_tmp_oppo(DIM+1),vector_tmp_oppo(DIM+2),vector_tmp_oppo(1:H_hr.Dim),'sym');
else
fprintf('The homecell hamilton is not hermi, never mind, we will not hermi it enforcely!\n');
end
end
continue;
end
if j == 0
FORMAT = strcat("The opposite vector  hamilton ",fold(@strcat,repmat("
Input = num2cell(vector_tmp);
Input{DIM+3} = string(H_hr.HcoeL(i)');
fprintf(FORMAT, Input{:});
H_hr_tmp = H_hr_tmp.set_hop(H_hr.HcoeL(i)',vector_tmp_oppo(DIM+1),vector_tmp_oppo(DIM+2),vector_tmp_oppo(1:H_hr.Dim),'sym');
continue;
end
if ~isequal(H_hr.HcoeL(i) ,H_hr_tmp.HcoeL(j)')
if H_hr.HcoeL(i)== sym(0)
FORMAT = strcat("The testing vector hamilton ",fold(@strcat,repmat("
Input = num2cell(vector_tmp);
Input{DIM+3} = string(H_hr.HcoeL(j)');
fprintf(FORMAT,Input{:});
H_hr_tmp = H_hr_tmp.set_hop(H_hr.HcoeL(j)',vector_tmp(DIM+1),vector_tmp(DIM+2),vector_tmp(1:H_hr.Dim),'sym');
elseif H_hr.HcoeL(j)== sym(0)
FORMAT = strcat("The opposite vector hamilton ",fold(@strcat,repmat("
Input = num2cell(vector_tmp);
Input{DIM+3} = string(H_hr.HcoeL(j)');
fprintf(FORMAT,Input{:});
H_hr_tmp = H_hr_tmp.set_hop(H_hr.HcoeL(i)',vector_tmp_oppo(DIM+1),vector_tmp_oppo(DIM+2),vector_tmp_oppo(1:H_hr.Dim),'sym');
else
FORMAT = strcat("check on the
"[i:
'find in the
Input1 = num2cell(vector_tmp);
Input2 = num2cell(vector_tmp_oppo);
Input = [{i},Input1,Input2,{j}];
fprintf(FORMAT,Input{:});
%                                         fprintf(['check on the %3d th NRPT,\n' ...
tmpsym = (H_hr.HcoeL(i)+H_hr.HcoeL(j)')/2;
H_hr_tmp = H_hr_tmp.set_hop(tmpsym,vector_tmp(DIM+1),vector_tmp(DIM+2),vector_tmp(1:H_hr.Dim),'sym');
H_hr_tmp = H_hr_tmp.set_hop(tmpsym',vector_tmp_oppo(DIM+1),vector_tmp_oppo(DIM+2),vector_tmp_oppo(1:H_hr.Dim),'sym');
end
end
end
otherwise
end
case 'num'
fprintf('For num Hr, The requirement is less strick!');
for i = 1:H_hr.NRPTS
vector_tmp = H_hr.vectorL(i,:);
vector_tmp_oppo = -vector_tmp;
[~,j]=ismember(vector_tmp_oppo,H_hr_tmp.vectorL,'rows');
fprintf('check %d th : NRPT,\n the vector is %d %d %d,\n the opposite vector is %d,%d,%d the %d th NRPT\n',i,...
vector_tmp(1),vector_tmp(2),vector_tmp(3),...
vector_tmp_oppo(1),vector_tmp_oppo(2),vector_tmp_oppo(3),j);
if i == j
if ~isequal(H_hr.HnumL(:,:,i) ,H_hr_tmp.HnumL(:,:,j)')
fprintf('The homecell hamilton is not hermi, never mind, we will hermi it enforcely!\n');
disp(H_hr.HnumL(:,:,i));
H_hr_tmp = H_hr_tmp.set_hop_mat(H_hr.HnumL(:,:,i)'/2+H_hr_tmp.HnumL(:,:,j)/2,[0,0,0],'set');
fprintf('change into\n');
disp(H_hr_tmp.HnumL(:,:,i));
end
continue;
end
if j == 0
fprintf('The opposite vector hamilton does not exist, build it!\n');
H_hr_tmp = H_hr_tmp.set_hop_mat(H_hr.HnumL(:,:,i)',-vector_tmp,'set');
disp(H_hr_tmp.HnumL(:,:,i)');
continue;
elseif ~isequal(H_hr.HnumL(:,:,i) ,H_hr_tmp.HnumL(:,:,j)')
fprintf('The opposite vector hamilton is not hermi, replace it by strong mat! \n');
N1 = nnz(H_hr.HnumL(:,:,i));
N2 = nnz(H_hr_tmp.HnumL(:,:,j)');
if N1 >= N2
fprintf('The %d th NRPT is stonger!\n',i);
disp(H_hr.HnumL(:,:,i));
H_hr_tmp = H_hr_tmp.set_hop_mat(H_hr.HnumL(:,:,i)',-vector_tmp,'set');
elseif N1 < N2
fprintf('The %d th NRPT is stonger!\n',j);
disp(H_hr_tmp.HnumL(:,:,j)');
H_hr_tmp = H_hr_tmp.set_hop_mat(H_hr_tmp.HnumL(:,:,j)',vector_tmp,'set');
else
end
elseif isequal(H_hr.HnumL(:,:,i),H_hr_tmp.HnumL(:,:,j)')
disp('hermi test pasts');
else
warning('!!!!!');
end
end
end
if giveback
H_hr = H_hr_tmp.rewind();
else
H_hr = H_hr_tmp;
end
end
