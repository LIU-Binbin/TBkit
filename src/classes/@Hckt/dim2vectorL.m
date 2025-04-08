function vectorL = dim2vectorL(dim,maxR)
if dim < 4 && maxR <= 1
switch dim
case 1
vectorL= [0;1];
case 2
vectorL= [0,0;1,0;0,1;1,1;1,-1];
case 3
vectorL= [...
0,0,0;...
1,0,0;...
0,1,0;...
0,0,1;...
1,1,0;...
0,1,1;...
1,0,1;...
-1,1,0;...
0,-1,1;...
1,0,-1;...
1,1,1;...
];
case 5
end
else
for i = 1:dim
VectorPre{i} = (-maxR:maxR).';
end
VectorStore = fold(@park.SemiProductVector,VectorPre);
VectorStore = sortrows(VectorStore,'descend');
AvectorL = zeros(1,dim);
BvectorL = AvectorL;
countA = 1;
countB = 1;
for i = 1:size(VectorStore,1)
vector = VectorStore(i,:);
if ismember(-vector,BvectorL,'row') || ismember(vector,BvectorL,'row')
continue
else
countA = countA + 1;
countB = countB + 1;
AvectorL(countA,:) = vector;
BvectorL(countB,:) = -vector;
end
end
vectorL = AvectorL;
end
end
