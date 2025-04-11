function HcktObj = HR2Hckt(H_hr,options,options_homecell,options_para)
% HR2HCKT Convert HR object to Hckt circuit representation
%
%   HcktObj = HR2HCKT(H_hr,options,options_homecell,options_para)
%   converts a tight-binding Hamiltonian to an equivalent circuit model.
%
%   INPUT ARGUMENTS:
%       H_hr - HR object to convert
%       options - Main conversion options:
%           title: Circuit title (default: 'HcktFromHR')
%           dim: Dimension (default: 2)
%           autominus: Auto-adjust signs (logical, default: false)
%           magnitude: Component magnitude ('p','n','u','m')
%           mode: Conversion mode ('real' or 'sigma')
%       options_homecell - Home cell options
%       options_para - Parameter options
%
%   OUTPUT ARGUMENTS:
%       HcktObj - Resulting circuit object
%
%   NOTES:
%       - Supports both real and Pauli matrix representations
%       - Handles different component magnitude scales
%
%   SEE ALSO:
%       Hckt, Subckt
%
%   AUTHOR:
%       [Your Name] ([Your Email])
%       [Creation Date]
arguments
H_hr
options.title = 'HcktFromHR';
options.dim = 2;
options.autominus = false;
options.magnitude = 'p';
options.mode {mustBeMember(options.mode,{'real','sigma'})}= 'real';
options_homecell.homecell = "normal";
options_homecell.defaultparameters = false;
options_para.prefix = 100;
end
VarC0=1;
Var2C0=2;
VarR0=20;
VarR0_2=0.5;
switch options.magnitude
case 'p'
Cmagnitude = 'p';
Lmagnitude = 'u';
Rmagnitude = 'k';
case 'n'
Cmagnitude = 'n';
Lmagnitude = 'u';
Rmagnitude = '';
case 'u'
Cmagnitude = 'p';
Lmagnitude = 'u';
Rmagnitude = 'k';
case 'm'
Cmagnitude = 'p';
Lmagnitude = 'u';
Rmagnitude = 'k';
end
TITLE = options.title;
Dim = options.dim;
HomeVector = zeros(1,Dim);
NBAND= double(H_hr.WAN_NUM);
switch options.mode
case 'real'
H_hr = H_hr.ForceTolist();
if options.autominus
warning off;
for i = 1:H_hr.WAN_NUM
BaseOnsite = sum(H_hr.HnumL(H_hr.vectorL(:,H_hr.Dim+1)==i));
H_hr = H_hr.set_hop(BaseOnsite,...
i,i,zeros(1,H_hr.Dim),'set');
end
end
HcktObj = Hckt('title',TITLE,'Nports',NBAND,'vectorL',HomeVector,'magnitude',options.magnitude);
fprintf('Limitations: Only Numerical HR, real hopping support.\n');
HomeCell = Subckt.FromHomecellList(H_hr,'magnitude',options.magnitude,'prefix',options_para.prefix);
HcktObj = HcktObj.set_home(HomeCell,1:round(NBAND/2),round(NBAND/2)+1:H_hr.WAN_NUM);
Cplus  = Subckt('Xplus 1 2  C') ;
Cminus  = Subckt('Xminus 1 2  minusC') ;
vectorList = double(H_hr.vectorL);
for n = 1:H_hr.NRPTS
Rvector = vectorList(n,1:Dim);
if isequal(Rvector,HomeVector)
continue;
end
i = vectorList(n,H_hr.Dim+1);
j = vectorList(n,H_hr.Dim+2);
if H_hr.HnumL(n) > 0
HoppingSckt = Cminus;
elseif H_hr.HnumL(n) < 0
HoppingSckt = Cplus;
end
HcktObj = HcktObj.set_hop(Rvector,HoppingSckt,i,j,['C_hopping = ',num2str(abs(H_hr.HnumL(n)*options_para.prefix )),options.magnitude]);
end
case 'sigma'
H_hr = H_hr.ForceToMat();
BasisC3_origin   =Subckt('XBasisC3_origin   1 2 3 TOGND BasisC3_origin  '            ,'magicnumber',6);
HoppingDist{1}   =Subckt('XPlusSigma0       l1 l2 l3 r1 r2 r3 TOGND PlusSigma0      ','magicnumber',9);
HoppingDist{2}   =Subckt('XMinusSigma0      l1 l2 l3 r1 r2 r3 TOGND MinusSigma0     ','magicnumber',9);
HoppingDist{3}   =Subckt('XPlusiSigma0      l1 l2 l3 r1 r2 r3 TOGND PlusiSigma0     ','magicnumber',9);
HoppingDist{4}   =Subckt('XMinusiSigma0     l1 l2 l3 r1 r2 r3 TOGND MinusiSigma0    ','magicnumber',9);
HoppingDist{5}   =Subckt('XPlusSigma1       l1 l2 l3 r1 r2 r3 TOGND PlusSigma1      ','magicnumber',9);
HoppingDist{6}   =Subckt('XMinusSigam1      l1 l2 l3 r1 r2 r3 TOGND MinusSigma1     ','magicnumber',9);
HoppingDist{7}   =Subckt('XPlusiSigma1      l1 l2 l3 r1 r2 r3 TOGND PlusiSigma1     ','magicnumber',9);
HoppingDist{8}   =Subckt('XMinusiSigma1     l1 l2 l3 r1 r2 r3 TOGND MinusiSigma1    ','magicnumber',9);
HoppingDist{9}   =Subckt('XPlusGen3Sigma2   l1 l2 l3 r1 r2 r3 TOGND PlusGen3Sigma2  ','magicnumber',9);
HoppingDist{10}  =Subckt('XMinusGen3Sigma2  l1 l2 l3 r1 r2 r3 TOGND MinusGen3Sigma2 ','magicnumber',9);
HoppingDist{11}  =Subckt('XPlusiGen3Sigma2  l1 l2 l3 r1 r2 r3 TOGND PlusiGen3Sigma2 ','magicnumber',9);
HoppingDist{12}  =Subckt('XMinusiGen3Sigma2 l1 l2 l3 r1 r2 r3 TOGND MinusiGen3Sigma2','magicnumber',9);
HoppingDist{13}  =Subckt('XPlusGen3Sigma3   l1 l2 l3 r1 r2 r3 TOGND PlusGen3Sigma3  ','magicnumber',9);
HoppingDist{14}  =Subckt('XMinusGen3Sigma3  l1 l2 l3 r1 r2 r3 TOGND MinusGen3Sigma3 ','magicnumber',9);
HoppingDist{15}  =Subckt('XPlusiGen3Sigma3  l1 l2 l3 r1 r2 r3 TOGND PlusiGen3Sigma3 ','magicnumber',9);
HoppingDist{16}  =Subckt('XMinusiGen3Sigma3 l1 l2 l3 r1 r2 r3 TOGND MinusiGen3Sigma3','magicnumber',9);
if NBAND ~= 2
fprintf('Choosing wrong mode!\n');
error('Nband ~= 2!!!!!');
end
HcktObj = Hckt('title',TITLE,'Nports',NBAND+1,'vectorL',HomeVector,'magnitude',options.magnitude);
fprintf('Limitations: Only Numerical HR, Two band model support.\n');
switch options_homecell.homecell
case 'normal'
HomeCell = BasisC3_origin;
otherwise
HomeCell = BasisC3_origin;
end
OnsiteMat       = H_hr.HnumL(:,:,H_hr.Line_000);
OnsitePotencial = OnsiteMat(1,1);
parameters = ['VarCg = ',num2str(abs(OnsitePotencial*VarC0)),Cmagnitude];
ScktObjDevice = 'X';
ScktObjName = 'Pri';
ScktObjNode = [string(1:3),'TOGND'];
ScktObjNetlist = HomeCell.netlist;
ScktObjDescription = parameters ;
ScktObj = Subckt(ScktObjDevice,ScktObjName,ScktObjNode,ScktObjDescription,ScktObjNetlist);
HcktObj = HcktObj.set_home(ScktObj,1:3,4:6);
vectorList = double(H_hr.vectorL);
HnumList = double(H_hr.HnumL);
for n = 1:H_hr.NRPTS
Rvector = vectorList(n,1:Dim);
if isequal(Rvector,HomeVector)
continue;
end
[CoeForPauli] = TBkit.pauliDecompositionNumerial(HnumList(:,:,n));
reaLCoeForPauli = real(CoeForPauli);
imagCoeForPauli = imag(CoeForPauli);
Coe16 = zeros(1,16);
PlusReal  = [1,5,9,13];
MinusReal = [2,6,10,14];
PlusImag  = [3,7,11,15];
MinusImag = [4,8,12,16];
SelectL = reaLCoeForPauli<0;
if ~isempty(SelectL)
Coe16( PlusReal(SelectL)) =    reaLCoeForPauli(SelectL);
end
SelectL = reaLCoeForPauli>0;
if ~isempty(SelectL)
Coe16(MinusReal(SelectL>0)) =   -reaLCoeForPauli(SelectL);
end
SelectL = imagCoeForPauli<0;
if ~isempty(SelectL)
Coe16( PlusImag(SelectL)) =    imagCoeForPauli(SelectL);
end
SelectL = imagCoeForPauli>0;
if ~isempty(SelectL)
Coe16(MinusImag(SelectL)) =   -imagCoeForPauli(SelectL);
end
for i = 1:16
if Coe16(i) == 0
else
if i < 9
Gen3Factor = 1;
else
Gen3Factor = 1;
end
if options_homecell.defaultparameters
Parameters = [];
else
Parameters =  [
' VarC0 = ',num2str(abs(Coe16(i)*Gen3Factor*VarC0)),Cmagnitude,...
' Var2C0 = ',num2str(abs(Coe16(i)*Gen3Factor*Var2C0)),Cmagnitude,...
' VarR0 = ',num2str(abs(Coe16(i)*Gen3Factor*VarR0)),Rmagnitude,...
' VarR0_2 = ',num2str(abs(Coe16(i)*Gen3Factor*VarR0_2)),Rmagnitude ...
];
end
HcktObj = HcktObj.set_hop(Rvector,HoppingDist{i},[1 2 3],[1 2 3],Parameters);
end
end
end
case 'Gell_Mann'
end
end
