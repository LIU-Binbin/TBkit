function cutlist= unique_label2cutlist(unique_label,NRPTS)
cutlist(:,1)= unique_label;
cutlist(1:end-1,2)= unique_label(2:end)-1;
cutlist(end,2) = NRPTS;
end
