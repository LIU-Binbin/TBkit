function Portlist_gen(fid,AllPortStrL,options,options2)
arguments
fid = 1;
AllPortStrL = [""];
options.type = 'v';
options2.comment = false;
options2.probenode = 'Allnode';
options2.prefix = "+ ";
options2.father = "";
options2.nodelist = [];
end
if options2.comment
switch options2.probenode
case 'Allnode'
for i = 1:numel(AllPortStrL)
fprintf(fid,"*"+options2.prefix+options.type+"(%s)\n",AllPortStrL(i));
end
case 'Selectnode'
for i = 1:numel(AllPortStrL)
for j = 1:numel(options2.nodelist)
fprintf(fid,"*"+options2.prefix+options.type+"(%s"+options2.father+"%s)\n",AllPortStrL(i),options2.nodelist(j));
end
end
end
else
switch options2.probenode
case 'Allnode'
for i = 1:numel(AllPortStrL)
fprintf(fid,options2.prefix+options.type+"(%s)\n",AllPortStrL(i));
end
case 'Selectnode'
for i = 1:numel(AllPortStrL)
for j = 1:numel(options2.nodelist)
fprintf(fid,options2.prefix+options.type+"(%s"+options2.father+"%s)\n",AllPortStrL(i),options2.nodelist(j));
end
end
end
end
end
