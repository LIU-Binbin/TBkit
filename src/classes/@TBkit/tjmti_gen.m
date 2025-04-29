function TBkitObj = tjmti_gen(TBkitObj,mode)
%TJMTI_GEN Generate tj-ti terms for tight-binding model
%
%   Syntax:
%       TBkitObj = tjmti_gen(TBkitObj,mode)
%
%   Description:
%       Generates tj-ti terms and their Fourier transforms for the
%       tight-binding Hamiltonian. Can operate in numeric or symbolic mode.
%
%   Inputs:
%       TBkitObj - TBkit object
%       mode     - Operation mode ('num' or 'sym', default='num')
%
%   Output:
%       TBkitObj - Modified TBkit object with tjmti terms
if nargin < 2
    mode = 'num';
end
TBkitObj = TBkitObj.timtj_gen(mode);
% Prepare tj - ti
for i = 1:2
    TBkitObj.tjmti{i} = - TBkitObj.timtj{i};
end
if strcmp(mode,'sym')
    ExpInnerTerm = matrixtimespage(TBkitObj.VarsSeqLcart(1:TBkitObj.Dim),TBkitObj.tjmti{1});
    TBkitObj.tjmti{3} = exp(1i*(sum(ExpInnerTerm,3)));
    ExpInnerTermFrac = matrixtimespage(TBkitObj.VarsSeqLcart(1:TBkitObj.Dim),TBkitObj.tjmti{2});
    TBkitObj.tjmti{4} =  exp(1i*(sum(ExpInnerTermFrac,3)));
else
    TBkitObj.tjmti{3} = [];
    TBkitObj.tjmti{4} = [];
end
end