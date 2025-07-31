function [klist_dos,klist1,klist2]=kmesh3D(TBkitObj,mesh,kz,mode,options)
%Kmesh3D Generate 3D k-mesh for DOS, band structure, or Fermi arc calculations.
%   [klist_dos, klist1, klist2] = Kmesh3D(TBkitObj, mesh, kz, mode, options)
%   generates k-point grids in reciprocal space based on specified modes.
%
%   This function supports three modes: 'bulk-DOS' for 3D uniform grids,
%   '3Dbandplot' for 2D slices at fixed kz, and 'fermiarc' for Fermi arc paths.
%
% Inputs:
%   TBkitObj   - Object containing structural information (Rm and Gk matrices).
%   mesh       - Grid dimensions [Nx, Ny] or [Nx, Ny, Nz] (default: [21,21]).
%   kz         - Scalar or vector for kz values (required in '3Dbandplot' mode).
%   mode       - Operation mode: 'bulk-DOS', '3Dbandplot', or 'fermiarc' (default: 'bulk-DOS').
%   options    - Optional structure with fields:
%       KPOINTS - Filename to export k-points (default: '').
%
% Outputs:
%   klist_dos  - Generated k-points matrix (N x 3).
%   klist1     - First dimension k-grid coordinates (varies by mode).
%   klist2     - Second dimension k-grid coordinates (varies by mode).
%
% Examples:
%   % Generate bulk DOS grid
%   [klist] = kmesh3D(TB, [30,30,30]);
%
%   % Generate 2D band slice at kz=0.5
%   [klist, kx, ky] = kmesh3D(TB, [50,50], 0.5, '3Dbandplot');
%
%   % Generate Fermi arc points and export to file
%   opts.KPOINTS = 'kpoints.dat';
%   kmesh3D(TB, [100,100], kfermi, 'fermiarc', opts);
arguments
    TBkitObj ;
    mesh = [21,21];
    kz = [];
    mode = 'bulk-DOS';
    options.KPOINTS = '';
end
%--------  nargin  --------
if isempty(kz) && strcmp(mode,'bulk-DOS')
elseif strcmp(mode,'bulk-DOS')
    mode = '3Dbandplot';
end
%--------  init  --------
%             TBkitObj.Gk  = (eye(3)*2*pi/TBkitObj.Rm)' ;% Reciprocal vector
%--------  distri  --------
switch mode
    case 'bulk-DOS'
        b_x=sum(TBkitObj.Gk(:,1));
        b_y=sum(TBkitObj.Gk(:,2));
        b_z=sum(TBkitObj.Gk(:,3));
        xnodes=mesh(1);
        ynodes=mesh(2);
        znodes=mesh(3);
        xlist=linspace(-b_x,b_x,xnodes);
        ylist=linspace(-b_y,b_y,ynodes);
        zlist=linspace(-b_z,b_z,znodes);
        klist_dos=[];
        for i=1:xnodes
            for j=1:ynodes
                for k=1:znodes
                    klist_temp=[xlist(i) ylist(j) zlist(k)];
                    klist_dos=[klist_dos;klist_temp];
                end
            end
        end
        klist1 = xlist;
        klist2 = ylist;
    case '3Dbandplot'
        disp('Primitive cell')
        disp(TBkitObj.Rm);
        disp('Reciprocal Lattice')
        disp(TBkitObj.Gk);
        disp("Rm*Gk' ");
        disp(sym(TBkitObj.Rm*TBkitObj.Gk'));
        xnodes=mesh(1);
        ynodes=mesh(2);
        klist1=linspace(-0.5,0.5,xnodes);
        klist2=linspace(-0.5,0.5,ynodes);
        %     kz_r = kz*Gk(3,:);
        klist_dos=[];
        for i=1:xnodes
            for j=1:ynodes
                klist_temp=[klist1(i) klist2(j) kz]*TBkitObj.Gk;
                klist_dos=[klist_dos;klist_temp];
            end
        end
    case 'fermiarc'
        %disp('fermi-arc');
        xnodes = mesh(1);
        ynodes = mesh(2);
        kfermiarc = kz;
        klist1 = [linspace(0,kfermiarc(2,1),xnodes)',...
            linspace(0,kfermiarc(2,2),xnodes)',...
            linspace(0,kfermiarc(2,3),xnodes)',...
            ];
        klist2 = [linspace(0,kfermiarc(3,1),ynodes)',...
            linspace(0,kfermiarc(3,2),ynodes)',...
            linspace(0,kfermiarc(3,3),ynodes)',...
            ];
        klist_1 = klist1 ;
        klist_2 = klist2 ;
        %         disp( klist_1 );
        %         disp( klist_2 );
        klist_dos = zeros(xnodes*ynodes,3);
        count =0;
        for i=1:xnodes
            for j=1:ynodes

                count = count+1;
                klist_dos(count,:) = klist_1(i,:)+klist_2(j,:)+ kfermiarc(1,:);
            end
        end
        klist1  = klist_1 + kfermiarc(1,:);
        klist2  = klist_2 + kfermiarc(1,:);
        %         klist1 = linspace(kfermiarc(1,1),max(kfermiarc(2,:)),xnodes)';
        %         klist2 = linspace(kfermiarc(1,2),max(kfermiarc(3,:)),xnodes)';
end
if strcmp(options.KPOINTS,'')
else
    fileID = fopen(options.KPOINTS,'w');
    for i = 1:size(klist1,1)
        fprintf(fileID,'%7.5f %7.5f %7.5f\n',klist2(1,1)+klist1(i,1),klist2(1,2)+klist1(i,2),klist2(1,3)+klist1(i,3));
        fprintf(fileID,'%7.5f %7.5f %7.5f\n',klist2(end,1)+klist1(i,1),klist2(end,2)+klist1(i,2),klist2(end,3)+klist1(i,3));
        fprintf(fileID,'\n');
    end
    fclose(fileID);
end
end

       