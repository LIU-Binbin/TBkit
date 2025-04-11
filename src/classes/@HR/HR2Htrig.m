function H_htrig = HR2Htrig(H_hr,options)

% HR2HTRIG Convert HR object to Htrig (trigonometric) representation
%
%   H_htrig = HR2HTRIG(H_hr,options) converts a real-space Hamiltonian
%   to trigonometric (Bloch) representation.
%
%   INPUT ARGUMENTS:
%       H_hr - HR object to convert
%       options - Conversion options:
%           sym: Use symbolic math (logical, default: true)
%           num: Force numerical conversion (logical, default: false)
%           fast: Use fast conversion mode (logical, default: false)
%           direct: Direct conversion (logical, default: false)
%           Type: Output type ('exp','sincos','mat','list')
%
%   OUTPUT ARGUMENTS:
%       H_htrig - Resulting Htrig object
%
%   NOTES:
%       - Supports multiple representation formats
%       - Includes progress reporting for large systems
%
%   SEE ALSO:
%       Htrig, HR2HK
%
%   AUTHOR:
%       [Your Name] ([Your Email])
%       [Creation Date]

arguments
H_hr HR;
options.sym =true;
options.num = false;
options.fast = false
options.direct = false;
options.Accuracy = 1e-6;
options.Type {mustBeMember(options.Type,{'exp','sincos','mat','list'})}='sincos' ;
end
if H_hr.num
options.num = true;
options.sym = false;
end
if  options.sym && ~options.fast
H_hr.Rm = sym(H_hr.Rm);
H_hr.orbL = sym(H_hr.orbL);
fprintf('please check whether the Rm and orbL are properly set.\n');
disp(H_hr.Rm);
disp(H_hr.orbL);
end
if options.direct
Hsym = sym(H_hr);
H_htrig = Htrig(Hsym,'Dim',H_hr.Dim);
H_htrig = H_htrig.TBkitCopy(H_hr);
elseif options.fast
if ~strcmp(H_hr.Type,'list')
H_hr = H_hr.rewrite('Accuracy',options.Accuracy);
else
H_hr = H_hr.simplify(options.Accuracy);
end
H_htrig = Htrig(H_hr.WAN_NUM,'Type','list','Dim',H_hr.Dim);
H_htrig = H_htrig.TBkitCopy(H_hr);
tiL = H_htrig.orbL(H_hr.vectorL(:,H_hr.Dim+1),:);
tjL = H_htrig.orbL(H_hr.vectorL(:,H_hr.Dim+2),:);
RtjmtiL = (double(H_hr.vectorL(:,1:H_hr.Dim))+ tjL-tiL);
ijL = double(H_hr.vectorL(:,[H_hr.Dim+1,H_hr.Dim+2]));
if options.num
H_htrig.HnumL = H_hr.HnumL;
H_htrig.HsymL_numL = [RtjmtiL*double(H_htrig.Rm),ijL] ;
H_htrig.num = true;
H_htrig.coe = false;
end
if options.sym
H_htrig.HcoeL = H_hr.HcoeL;
H_htrig.HsymL_coeL = [RtjmtiL*H_htrig.Rm,ijL] ;
H_htrig.coe = true;
H_htrig.num = false;
end
if ~strcmp(options.Type,'list')
H_htrig = H_htrig.rewrite();
end
else
H_htrig = Htrig(H_hr.WAN_NUM,'Type',options.Type,'Dim',H_hr.Dim );
H_htrig = H_htrig.TBkitCopy(H_hr);
if options.num
H_hr.HcoeL = sym(H_hr.HnumL);
end
H_htrig = H_htrig.tjmti_gen('sym');
VarsUsing = H_htrig.VarsSeqLcart(1:H_htrig.Dim);
vectorList = double(H_hr.vectorL);
pb = TBkit_tool_outer.CmdLineProgressBar('Transforming ');
if strcmp(H_hr.Type,'list')
tji_mat_cart = H_htrig.tjmti{1};
NRPTS_ = H_hr.NRPTS;
for k = 1:NRPTS_
pb.print(k,NRPTS_,' th hopping into Htrig obj ...');
SymHopping = H_hr.HcoeL(k);
if SymHopping ~= sym(0)
i = vectorList(k,H_hr.Dim+1);
j = vectorList(k,H_hr.Dim+2);
SymVar = exp(...
1i...
*VarsUsing*((...
vectorList(k,1:H_hr.Dim)*H_hr.Rm)'+reshape(tji_mat_cart(i,j,:),[3,1])...
));
H_htrig = H_htrig.set_hop(SymHopping,SymVar,[i,j]);
end
end
else
tji_mat_kcart = H_htrig.tjmti{3};
NRPTS_ = H_hr.NRPTS;
exp_preL = exp(1i...
*VarsUsing*((...
vectorList*H_hr.Rm)'));
WAN_WAN = numel(tji_mat_kcart);
sizeWAN = size(tji_mat_kcart);
sym0 = sym(0);
[iL,jL] = ind2sub(sizeWAN,1:WAN_WAN);
for i = 1:NRPTS_
pb.print(i,NRPTS_,' th H(NRPT) into Htrig obj ...');
SymHoppingMat = H_hr.HcoeL(:,:,i);
exp_preL_i  = exp_preL(i);
for j = 1:WAN_WAN
SymHopping = SymHoppingMat(j);
if SymHopping ~= sym0
exp_pre = exp_preL_i*tji_mat_kcart(j);
SymVar = simplify(rewrite(exp_pre,'sincos'));
H_htrig = H_htrig.set_hop(SymHopping,SymVar,[iL(j),jL(j)]);
end
end
end
end
pb.delete();
end
end
