function H_hr = pythtb_gen(H_hr,filename,options)
arguments
H_hr HR ;
filename char='HRTB.py'
options.sym = true;
options.band logical= true;
options.python_env = '/usr/bin/python';
options.KPOINTS = 'KPOINTS';
options.WilsonLoop = false;
end
H_hr = H_hr.rewrite();
fileID = fopen(filename,'w');
fprintf(fileID,['#!',options.python_env ,'\n']);
fprintf(fileID,'from pythtb import *\n');
fprintf(fileID,'# lattice vectors and orbital positions\n');
latstring = "lat=[";
for i = H_hr.Dim
latstring = latstring+TBkit.mat2str_python(H_hr.Rm(i,:))+", ";
end
latstring = latstring + "]";
fprintf(fileID,latstring+'\n');
orbstring = "orb=[\n";
for i =1:H_hr.WAN_NUM
orbstring = orbstring + TBkit.mat2str_python(H_hr.orbL(i,:))+", \n";
end
orbstring = orbstring + "]";
fprintf(fileID,orbstring+'\n');
fprintf(fileID,'\n');
fprintf(fileID,'# three-dimensional tight-binding model from HR\n');
fprintf(fileID,['HRTB=tb_model(',num2str(3),',',num2str(3),', lat, orb)\n']);
fprintf(fileID,'\n');
fprintf(fileID,'# define hopping between orbitals\n');
SelectLabel_home = all([H_hr.vectorL(:,1:H_hr.Dim) == 0,H_hr.vectorL(:,H_hr.Dim+2)>=H_hr.vectorL(:,H_hr.Dim+1)],2);
SelectLabel_other = all(H_hr.vectorL(:,1:H_hr.Dim) >=0,2) & ~all(H_hr.vectorL(:,1:H_hr.Dim) == 0,2);
SelectLabel = SelectLabel_home|SelectLabel_other;
SelectVector=H_hr.vectorL(SelectLabel,:);
if options.sym
SelectHop=H_hr.HcoeL(SelectLabel);
for i = 1:length(H_hr.symvar_list)
fprintf(fileID,string(SelectHop(i))+"= 1 \n");
end
else
SelectHop=H_hr.HnumL(SelectLabel);
end
for i = 1:numel(SelectHop)
fprintf(fileID,'HRTB.set_hop('+string(SelectHop(i))+","+...
num2str(SelectVector(i,H_hr.Dim+1)-1)+","+num2str(SelectVector(i,H_hr.Dim+2)-1)+","+...
TBkit.mat2str_python(SelectVector(i,1:H_hr.Dim))+')\n');
end
fprintf(fileID,'\n');
if options.band
[kpoints,nodes,kpoints_name_tmp] = TBkit.KPOINTS_read(options.KPOINTS);
nkline = length(kpoints_name_tmp)-1;
kpoints_f = kpoints(1+(0:nkline-1)*2,:);
kpoints_f = [kpoints_f;kpoints(end,:)];
fprintf(fileID,'# solve model on a path in k-space\n');
fprintf(fileID,'k=[');
for i = 1:size(kpoints_f,1)
fprintf(fileID,[TBkit.mat2str_python(kpoints_f(i,:)),',']);
end
fprintf(fileID,']\n');
fprintf(fileID,['(k_vec,k_dist,k_node)=HRTB.k_path(k,',num2str(nodes),')\n']);
fprintf(fileID,'evals=HRTB.solve_all(k_vec)\n');
fprintf(fileID,'# plot bandstructure\n');
fprintf(fileID,'import matplotlib.pyplot as plt\n');
fprintf(fileID,'fig, ax = plt.subplots()\n');
fprintf(fileID,['for i in range(1,',num2str(H_hr.WAN_NUM+1),'):\n']);
fprintf(fileID,'    ax.plot(k_dist,evals[i-1,:])\n');
fprintf(fileID,'ax.set_xticks(k_node)\n');
fprintf(fileID,['ax.set_xticklabels(',TBkit.mat2str_python(kpoints_name_tmp),')\n']);
fprintf(fileID,'ax.set_xlim(k_node[0],k_node[-1])\n');
fprintf(fileID,'fig.savefig("band.png")\n');
end
if options.WilsonLoop
fprintf(fileID,'#calculate my-array\n');
fprintf(fileID,'my_array=wf_array(HRTB,[41,41,41])\n');
fprintf(fileID,'# solve model on a regular grid, and put origin of\n');
fprintf(fileID,'# Brillouin zone at [0,0,0]  point\n');
fprintf(fileID,'my_array.solve_on_grid([0,0,0])\n');
fprintf(fileID,'\n');
fprintf(fileID,'# calculate Berry phases around the BZ in the k_x direction\n');
fprintf(fileID,'# (which can be interpreted as the 1D hybrid Wannier centers\n');
fprintf(fileID,'# in the x direction) and plot results as a function of k_y\n');
fprintf(fileID,'#\n');
fprintf(fileID,'# Following the ideas in\n');
fprintf(fileID,'#   A.A. Soluyanov and D. Vanderbilt, PRB 83, 235401 (2011)\n');
fprintf(fileID,'#   R. Yu, X.L. Qi, A. Bernevig, Z. Fang and X. Dai, PRB 84, 075119 (2011)\n');
fprintf(fileID,'# the connectivity of these curves determines the Z2 index\n');
fprintf(fileID,'#\n');
fprintf(fileID,['wan_cent = my_array.berry_phase(',...
TBkit.mat2str_python(1:H_hr.WAN_NUM/2),...
',dir=1,contin=False,berry_evals=True)\n']);
fprintf(fileID,'wan_cent/=(2.0*np.pi)\n');
fprintf(fileID,'\n');
fprintf(fileID,'nky=wan_cent.shape[0]\n');
fprintf(fileID,'ky=np.linspace(0.,1.,nky)\n');
fprintf(fileID,'# draw Wannier center positions\n');
fprintf(fileID,'fig, ax2 = plt.subplots()\n');
fprintf(fileID,['for iband in range(0,',num2str(H_hr.WAN_NUM/2),'):\n']);
fprintf(fileID,'    ax2.plot(ky,wan_cent[iband,:,0],"k.")\n');
fprintf(fileID,'ax2.set_ylim(-1.0,1.0)\n');
fprintf(fileID,'ax2.set_ylabel(''Wannier center along x'')\n');
fprintf(fileID,'ax2.set_xlabel(''k_y'')\n');
fprintf(fileID,'ax2.set_xticks([0.0,0.5,1.0])\n');
fprintf(fileID,'ax2.set_xlim(0.0,1.0)\n');
fprintf(fileID,'ax2.set_xticklabels([r"$0$",r"$\\pi$", r"$2\\pi$"])\n');
fprintf(fileID,'ax2.axvline(x=.5,linewidth=0.5, color=''k'')\n');
fprintf(fileID,'ax2.set_title("1D Wannier centers phase")\n');
fprintf(fileID,'\n');
fprintf(fileID,'fig.tight_layout()\n');
fprintf(fileID,'fig.savefig("wcc.pdf")\n');
end
end
