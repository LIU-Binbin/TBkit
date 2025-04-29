function fid = kpath_card_gen(TBkitobj,mode,filename)
%KPATH_CARD_GEN Generate k-path input file for WT code
%   fid = kpath_card_gen(TBkitobj, mode, fname) creates WT-compatible input
%
%   Inputs:
%       TBkitobj - Object containing k-path data
%       mode     - File mode ('wt' for WT code)
%       filename - Output filename (default: 'wt.in')
%
%   Output:
%       fid - File identifier
%
%   Currently only supports 'wt' mode
arguments
    TBkitobj ;
    mode = 'wt';
    filename = 'wt.in';
end
% function string
if strcmp(mode,'wt')
    fid = fopen(filename,'a');
    head_string="KPATH_BULK            ! k point path ";
    fprintf(fid,"%s\n",head_string);
    fprintf(fid,"%d\n",abs(length(TBkitobj.kpoints_name-1)));
    for i=1:abs(length(TBkitobj.kpoints_name-1))
        fprintf(fid,"%1s ",strrep(TBkitobj.kpoints_name(i),'GAMMA','G'));
        fprintf(fid,"%9f ",num2str(TBkitobj.kpoints_frac(2*i-1,1)));
        fprintf(fid,"%9f ",num2str(TBkitobj.kpoints_frac(2*i-1,2)));
        fprintf(fid,"%9f ",num2str(TBkitobj.kpoints_frac(2*i-1,3)));
        fprintf(fid,"%1s ",strrep(TBkitobj.kpoints_name(i+1),'GAMMA','G'));
        fprintf(fid,"%9s ",num2str(TBkitobj.kpoints_frac(2*i,1)));
        fprintf(fid,"%9s ",num2str(TBkitobj.kpoints_frac(2*i,2)));
        fprintf(fid,"%9s \n",num2str(TBkitobj.kpoints_frac(2*i,3)));
    end
    fclose(fid);
end
end
