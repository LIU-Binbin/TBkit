function M1 = Mirror_construct(options)
% MIRROR_CONSTRUCT Construct mirror symmetry operator (Mx, My, or Mz)
%
%   M1 = Mirror_construct(options)
%
%   DESCRIPTION:
%       Constructs the mirror symmetry operator Mx, My, or Mz, acting on 
%       a tight-binding basis set (position + orbital + spin) using 
%       data from a wannier90.wout file or an HR object.
%
%   INPUT: 
%       options.mode           - 'spinful' (default) or 'spinless'
%       options.Mirror         - 'Mx', 'My', or 'Mz' (default: 'Mz')
%       options.source         - 'wout' (default) or 'HR'
%       options.Mirror_plane   - 3-element vector, mirror plane coordinates 
%                                (default: [0, 0, 0.5])
%       options.filename       - path to .wout file (default: 'wannier90.wout')
%       options.POSCAR_name    - path to POSCAR (default: 'POSCAR')
%       options.HR             - HR object if source = 'HR'
%
%   OUTPUT:
%       M1 - Oper object, with M1.U as the mirror operator matrix
%
%   FEATURES:
%       - Handles orbital parity (from orbital_matrix)
%       - Supports spinless or spinful basis (including -i*sigma mirror ops)
%       - Auto-detects spin from quantumL column 4
%
%   SEE ALSO:
%       wout_read, HR class, Oper class


arguments
    options.mode = 'spinful'
    options.Mirror  = "Mz"
    options.source  = "wout"
    options.Mirror_plane double = [0,0,0.5]
    options.filename = 'wannier90.wout'
    options.POSCAR_name = 'POSCAR'
    %options.UsePOSCAR_Cordinates = true
    options.HR = HR
end

% ------------ Data loading -------------
if options.source == "wout"
    %opts_sub.UsePOSCAR_Cordinates = options.UsePOSCAR_Cordinates;
    [orbL, elementL, quantumL] = wout_read(options.filename, options.POSCAR_name);
else % HR source
    if isempty(options.HR)
        error('HR object must be provided when source = "HR"');
    end
    HR1=options.HR;
    orbL     = HR1.orbL;
    elementL = HR1.elementL;
    quantumL = HR1.quantumL;
end

%----------检查basis是否有spin------
if strcmp(options.mode,'spinful')
    if size(quantumL,2)<4 || abs(quantumL(1,4)) ~= 0.5
        fprintf('未检测quantumL基的自旋成分，已切换为spinless模式')
        options.mode='spinless';
    end
end

% --------------orbit table of Mirror----------------
% 数字矩阵，列分别为：
% [ l, m, Mx, My, Mz ]
orbital_matrix = [
     0,  0, +1, +1, +1;
     1,  0, +1, +1, -1;
     1, +1, -1, +1, +1;
     1, -1, +1, -1, +1;
     2,  0, +1, +1, +1;
     2, +1, -1, +1, -1;
     2, -1, +1, -1, -1;
     2, +2, +1, +1, +1;
     2, -2, -1, -1, +1;
     3,  0, +1, +1, -1;
     3, +1, -1, +1, -1;
     3, -1, +1, -1, -1;
     3, -3, -1, -1, -1;
     3, +3, -1, +1, +1;
     3, -2, +1, -1, -1;
     3, +2, -1, -1, -1;
];


matA = [orbL, elementL, quantumL];
matB = [];

Mx_spin = complex(0,-1)*[0 1;1 0];
My_spin = complex(0,-1)*[0 complex(0,-1);complex(0,1) 0];
Mz_spin = complex(0,-1)*[1 0;0 -1];
M1 = Oper();
M1.U  = zeros(size(orbL,1), size(orbL,1));
M_quantum = zeros(1, size(orbL,1));

%------构造Mirror映射后的矩阵
if strcmp(options.Mirror, 'Mz')
    M1.R  = diag([1,1,-1]);
    M1.Rf = diag([1,1,-1]);
    Mz_plane = options.Mirror_plane(3);
    for j=1:size(matA,1)
        line_tmp = matA(j,:);
        line_tmp(3) = 2 * Mz_plane - matA(j,3);
        matB = [matB; line_tmp];
    end
elseif strcmp(options.Mirror, 'Mx')
    M1.R  = diag([-1,1,1]);
    M1.Rf = diag([-1,1,1]);
    Mx_plane = options.Mirror_plane(1);
    for j=1:size(matA,1)
        line_tmp = matA(j,:);
        line_tmp(1) = 2 * Mx_plane - matA(j,1);
        matB = [matB; line_tmp];
    end
elseif  strcmp(options.Mirror, 'My')
    M1.R  = diag([1,-1,1]);
    M1.Rf = diag([1,-1,1]);
    My_plane = options.Mirror_plane(2);
    for j=1:size(matA,1)
        line_tmp = matA(j,:);
        line_tmp(2) = 2 * My_plane - matA(j,2);
        matB = [matB; line_tmp];
    end
else
    fprintf('mode should be choosed in Mx,My,Mz')
end

%-----------------------轨道部分的M本征值
for i = 1:size(orbL,1)
    for j = 1:size(orbital_matrix,1)
        if quantumL(i,2)== orbital_matrix(j,1) && quantumL(i,3)==orbital_matrix(j,2)
            if strcmp(options.Mirror, 'Mx')
                M_quantum(i)=orbital_matrix(j,3);
            elseif strcmp(options.Mirror, 'My')
                M_quantum(i)=orbital_matrix(j,4);
            elseif strcmp(options.Mirror, 'Mz')
                M_quantum(i)=orbital_matrix(j,5);
            end
            break
        end
    end
end

if strcmp(options.mode,'spinless')
    for i = 1:size(orbL,1)
            tol = 1e-6;
            rowA = matA(i,:);
            found = false;
            for ii = 1:size(matB,1)
                if all(abs(matB(ii,1:7) - rowA(1:7)) < tol)
                    j = ii;
                    found = true;
                    break;
                end
            end
            if ~found
                fprintf('error, some orbit not found their Mirror couple\n');
            end
            M1.U(i,j) = M_quantum(i);
    end

elseif strcmp(options.mode,'spinful') && strcmp(options.Mirror, 'Mz')
    for i = 1:size(orbL,1)
            tol = 1e-6;
            rowA = matA(i,:);
            found = false;
            for ii = 1:size(matB,1)  % 找相同自旋但是位置是M映射的
                if all(abs(matB(ii,:) - rowA) < tol)
                    j = ii;    
                    found = true;
                    break;
                end
            end
            if ~found
                fprintf('error, some orbit not found their Mirror couple\n');
            end
            if rowA(8) > 0
                M1.U(i,j) = M_quantum(i)*Mz_spin(1,1);  
            else
                M1.U(i,j) = M_quantum(i)*Mz_spin(2,2);
            end
    end

elseif strcmp(options.mode,'spinful') && strcmp(options.Mirror, 'Mx')
    for i = 1:size(orbL,1)
            tol = 1e-6;
            rowA = matA(i,:);
            found = false;
            for ii = 1:size(matB,1)  % 找相反自旋但是位置是M映射的
                if all(abs(matB(ii,1:7) - rowA(1:7)) < tol) && abs( matB(ii,8) + rowA(8) ) <tol
                    j = ii;    
                    found = true;
                    break;
                end
            end
            if ~found
                fprintf('error, some orbit not found their Mirror couple\n');
            end
            if rowA(8) > 0    %up to dn
                M1.U(i,j) = M_quantum(i)*Mx_spin(2,1);  
            else                  % dn to up
                M1.U(i,j) = M_quantum(i)*Mx_spin(1,2);
            end
    end

elseif strcmp(options.mode,'spinful') && strcmp(options.Mirror, 'My')
    for i = 1:size(orbL,1)
            tol = 1e-6;
            rowA = matA(i,:);
            found = false;
            for ii = 1:size(matB,1)  % 找相反自旋但是位置是M映射的
                if all(abs(matB(ii,1:7) - rowA(1:7)) < tol) && abs( matB(ii,8) + rowA(8) ) <tol
                    j = ii;    
                    found = true;
                    break;
                end
            end
            if ~found
                fprintf('error, some orbit not found their Mirror couple\n');
            end
            if rowA(8) > 0    %up to dn
                M1.U(i,j) = M_quantum(i)*My_spin(2,1);  
            else                  % dn to up
                M1.U(i,j) = M_quantum(i)*My_spin(1,2);
            end
    end

end
    
        
end

