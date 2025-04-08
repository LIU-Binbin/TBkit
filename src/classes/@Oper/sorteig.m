function [Asort,Usort] = sorteig(U,A)
if nargin <2
mode = 'eigenval';
NUM_WAN = length(U);
else
mode = 'whole';
NUM_WAN = length(A);
end
NBANDS = length(U);
if strcmp(mode,'whole')
if isa(U,'sym')
Asort=sym(zeros(NUM_WAN ,NBANDS ));
SortTmp=diag(simplify(U));
else
Asort=zeros(NUM_WAN ,NBANDS );
SortTmp=diag(U);
end
[Usort,IJ]=sort(double(SortTmp),1,'ComparisonMethod','real');
for jj=1:NBANDS
Asort(:,jj)=A(:,IJ(jj));
end
Usort = diag(Usort);
elseif strcmp(mode,'eigenval')
SortTmp=diag(U);
[Usort,~]=sort(SortTmp,1,'ComparisonMethod','real');
Usort = diag(Usort);
Asort = [];
end
end
