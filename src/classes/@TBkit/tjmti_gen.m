function TBkitObj = tjmti_gen(TBkitObj,mode)
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