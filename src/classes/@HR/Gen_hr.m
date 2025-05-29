function hrdat = Gen_hr(H_hr,filename,mode)
% GEN_HR Generate Wannier90 _hr.dat file from HR object
%
%   hrdat = GEN_hr(H_hr,filename,mode) writes Hamiltonian data to
%   Wannier90 format file.
%
%   INPUT ARGUMENTS:
%       H_hr - HR object to export
%       filename - Output file path (default: 'wannier90_hr.dat')
%       mode - Output mode (currently only 'hr_dat' supported)
%
%   OUTPUT ARGUMENTS:
%       hrdat - File handle of generated file
%
%   NOTES:
%       - Converts to dense matrix format if necessary
%       - Includes progress reporting
%       - Follows Wannier90 file format specification
%
%   SEE ALSO:
%       HR, full
%
%   AUTHOR:
%       [Your Name] ([Your Email])
%       [Creation Date]

if nargin < 3
    mode = 'hr_dat';
end
if nargin <2
    filename = 'wannier90_hr.dat';
end
if strcmp(H_hr.Type,'sparse')
    H_hr  = H_hr.full();
end
if strcmp(H_hr.Type,'list')
    H_hr  = H_hr.rewind();
end
hrdat=fopen(filename ,'w');
pb = CmdLineProgressBar('Writing hr.dat - NRPT:');
switch mode
    case 'hr_dat'
        date_=date;
        fprintf(hrdat,"%s\n",date_);
        fprintf(hrdat,"         %d\n",H_hr.WAN_NUM);
        fprintf(hrdat,"         %d\n",H_hr.NRPTS);
        COUNT_FLAG=0;
        for i =1:H_hr.NRPTS
            fprintf(hrdat,"    %d",1);
            COUNT_FLAG=COUNT_FLAG+1;
            if COUNT_FLAG == 15
                COUNT_FLAG = 0;
                fprintf(hrdat,"\n");
            end
        end
        fprintf(hrdat,"\n");
        for i=1:H_hr.NRPTS
            pb.print(i,H_hr.NRPTS,' ...');
            %fprintf('Wrinting (%d/%d)  NRPT\n',i,H_hr.NRPTS);
            for k=1:H_hr.WAN_NUM
                %fprintf('Wrinting %d th WAN_ORB \n',k);
                for j=1:H_hr.WAN_NUM
                    fprintf(hrdat,"    %d    %d    %d",H_hr.vectorL(i,1),H_hr.vectorL(i,2),H_hr.vectorL(i,3));
                    fprintf(hrdat,"    %d    %d",j,k);
                    real_part=real(H_hr.HnumL(j,k,i));
                    imag_part=imag(H_hr.HnumL(j,k,i));
                    fprintf(hrdat,"    %.16f    %.16f\n",real_part,imag_part);
                end
            end
        end
        fclose(hrdat);
    otherwise
end
end
