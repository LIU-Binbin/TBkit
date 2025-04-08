function H_hr = enlarge(H_hr,dir,amp)
if isa(dir,'char')
H_hr = H_hr.rewrite;
vectorList = H_hr.vectorL;
orbList= H_hr.orbL;
orbListdiff = orbList(vectorList(:,H_hr.Dim+1),:) - orbList(vectorList(:,H_hr.Dim+2),:);
orbListdiff_r = orbListdiff*H_hr.Rm;
switch dir
case 'x'
selectL = logical(orbListdiff_r(:,1));
case 'y'
selectL = logical(orbListdiff_r(:,2));
case 'z'
selectL = logical(orbListdiff_r(:,3));
end
H_hr.HnumL(selectL)=H_hr.HnumL(selectL)*amp;
else
switch size(dir,2)
case 1
selectL = H_hr.vectorL(:,dir) ;
case 2
case 3
case 5
end
end
end
