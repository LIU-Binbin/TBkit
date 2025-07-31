function [EIGENCAR_3D,klist1,klist2,WEIGHTCAR_3D,WAVECAR_3D] = EIGENCAR_gen_3D(H_hr,kmesh,k3d,options)
% EIGENCAR_GEN_3D Generate 3D eigenstate data from HR Hamiltonian
%
%   [EIGENCAR_3D,klist1,klist2,WEIGHTCAR_3D,WAVECAR_3D] = EIGENCAR_GEN_3D(H_hr,kmesh,k3d,options)
%   calculates eigenstates in 3D k-space for the given Hamiltonian.
%
%   INPUT ARGUMENTS:
%       H_hr - Hamiltonian in HR format
%       kmesh - Mesh dimensions [Nk1,Nk2] (default: [100,100])
%       k3d - 3D k-space basis vectors (4x3 matrix, default provided)
%       options - Structure with calculation options:
%           output: 'raw_data' or 'refined' (default: 'raw_data')
%           LWAVE: Calculate wavefunctions (logical, default: false)
%           cartisian: Use cartesian coordinates (logical, default: true)
%           fin_dir: Final direction for refined output (1-3, default: 3)
%           WEIGHTCAR: Calculate weight data (logical, default: false)
%           norb: Number of orbitals (default: -1)
%           fermi: Fermi level (default: 0)
%
%   OUTPUT ARGUMENTS:
%       EIGENCAR_3D - Eigenvalue data (3D array)
%       klist1,klist2 - k-point mesh coordinates
%       WEIGHTCAR_3D - Weight data (if requested)
%       WAVECAR_3D - Wavefunction data (if requested)
%
%   SEE ALSO:
%       HR, EIGENCAR_gen
%
%   AUTHOR:
%       [Your Name] ([Your Email])
%       [Creation Date]
arguments
    H_hr HR
    kmesh = [100,100];
    k3d = [-0.5, -0.5, 0; 1, 0, 0; 0, 1, 0];
    options.output ='raw_data';
    options.LWAVE = false;
    options.cartisian = true;
    options.fin_dir = 3;
    options.ProjectionMethod = 'slab';
    options.ProjectionStruct = struct('discrimination',0.1,'center',[0.5,0.5,0.5],'orientation',2,'sign',false);
    options.WEIGHTCAR = false;
    options.norb = -1;
    options.fermi = 0;
end
options2 = rmfield(options,{'output','cartisian','fin_dir'});
optionsCell = namedargs2cell(options2);
[~, H_hr.klist_frac, klist1, klist2] = kmeshgen(H_hr_arc.Rm, [k3d; 0 0 0], 'Nk1',kmesh(1), 'Nk2',kmesh(2));
if options.WEIGHTCAR
    [EIGENCAR_3D,WAVECAR_3D,WEIGHTCAR_3D] = H_hr.EIGENCAR_gen(...
        optionsCell{:});
else
    [EIGENCAR_3D,WAVECAR_3D] = H_hr.EIGENCAR_gen('LWAVE',options.LWAVE);
    WEIGHTCAR_3D = [];
end
if options.LWAVE
else
    WAVECAR_3D = [];
end
if options.cartisian
    options.output = 'refined';
    klist = H_hr.klist_frac*H_hr.Gk;
end
if options.LWAVE
    WEIGHTCAR_3D = WAVECAR_3D;
end
if strcmp(options.output ,'refined')
    switch options.fin_dir
        case 1
            klist1 =klist(:,2);
            klist2 =klist(:,3);
        case 2
            klist1 =klist(:,1);
            klist2 =klist(:,3);
        case 3
            klist1 =klist(:,1);
            klist2 =klist(:,2);
    end
    klist1 = reshape(klist1,kmesh);
    klist2 = reshape(klist2,kmesh);
    EIGENCAR_3D = reshape(EIGENCAR_3D.',kmesh(1),kmesh(2),[]);
    if options.WEIGHTCAR
        WEIGHTCAR_3D = reshape(WEIGHTCAR_3D.',kmesh(1),kmesh(2),[]);
    end
end
end
