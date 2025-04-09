function [vector_list,Coeffs_list] = Hstr_mapping(H_str,R_struct,mode)
%HSTR_MAPPING Map k-space terms to real-space hopping vectors/coefficients
%
% Syntax:
%   [v,c] = Hstr_mapping(H_str,R_struct)
%   [v,c] = Hstr_mapping(H_str,R_struct,mode)
%
% Inputs:
%   H_str - String representing k-space term (e.g., 'x', 'x*y')
%   R_struct - Struct with lattice constants (a,b,c fields)
%   mode - Coordinate system type:
%       'T' - Tetragonal (default)
%       'H' - Hexagonal
%
% Outputs:
%   vector_list - Nx3 array of hopping vectors
%   Coeffs_list - Corresponding coefficients for each vector
%
% Description:
%   Translates k-space polynomial terms (like k_x) to their real-space
%   tight-binding equivalents through Taylor expansion. Handles:
%   - First-order terms (x,y,z)
%   - Second-order terms (X,Y,Z)
%   - Cross terms (x*y)
%   - Composite terms through recursive decomposition
%
% Supported Terms:
%   Tetragonal ('T'):
%     1, x, y, z, X, Y, Z
%   Hexagonal ('H'):
%     Additional complex phase factors for hexagonal symmetry
%
% Example:
%   [v,c] = Hstr_mapping('x*y',struct('a',1,'b',1),'H')
if nargin <3
    mode = 'T';
end
if strcmp(mode,'T')
    switch H_str
        case '1'
            vector_list =  [0 0 0];
            Coeffs_list = 1;
        case 'x'
            vector_list =  [1 0 0;-1 0 0];
            Coeffs_list = (1/(2*R_struct.a*1i))*[1;-1];
        case 'y'
            vector_list =  [0 1 0;0 -1 0];
            Coeffs_list = (1/(2*R_struct.b*1i))*[1;-1];
        case 'z'
            vector_list =  [0 0 1;0 0 -1];
            Coeffs_list = (1/(2*R_struct.c*1i))*[1;-1];
        case 'X'
            vector_list =  [0 0 0;1 0 0;-1 0 0];
            Coeffs_list = (1/R_struct.a^2)*[2;-1;-1];
        case 'Y'
            vector_list =  [0 0 0;0 1 0;0 -1 0];
            Coeffs_list = (1/R_struct.b^2)*[2;-1;-1];
        case 'Z'
            vector_list =  [0 0 0;0 0 1;0 0 -1];
            Coeffs_list = (1/R_struct.c^2)*[2;-1;-1];
        otherwise
            disp(H_str);
            warning('the above Term cant be kp2tb, the result is wrong, please check!');
    end
elseif strcmp(mode,'H')
    switch H_str
        case '1'
            vector_list =  [0 0 0];
            Coeffs_list = 1;
        case 'x'
            vector_list =  [...
                1 0 0 ;-1  0 0;...
                0 1 0 ;0 -1 0;...
                -1 -1 0;1  1 0;...
                ];
            Coeffs_list = (1/(6*R_struct.a*1i))*[...
                2;-2;
                (exp(-pi*2i/3)+exp(pi*2i/3));-(exp(-pi*2i/3)+exp(pi*2i/3));
                (exp(-pi*4i/3)+exp(pi*4i/3));-(exp(-pi*4i/3)+exp(pi*4i/3));
                ];
        case 'y'
            vector_list =  [...
                0 1 0 ;0 -1 0;...
                -1 -1 0;1  1 0;...
                ];
            Coeffs_list = (1/(6*R_struct.b*1i))*[...
                -(exp(pi*1i/6)+exp(-pi*1i/6));(exp(pi*1i/6)+exp(-pi*1i/6));
                -(exp(pi*5i/6)+exp(-pi*5i/6));(exp(pi*5i/6)+exp(-pi*5i/6));
                ];
        case 'z'
            vector_list =  [0 0 1;0 0 -1];
            Coeffs_list = (1/2*R_struct.c*1i)*[1;-1];
        case 'X'
            vector_list =  [...
                0 0 0;1 0 0 ;-1  0 0;...
                0 0 0;1  1 0;-1 -1 0;...
                0 0 0;0 1 0 ;0 -1 0;...
                ];
            Coeffs_list = (1/(6*R_struct.a^2))*[...
                12;-6;-6;...
                0;0;0;...
                0;0;0;...
                ];
        case 'Y'
            vector_list =  [...
                0 0 0;1 0 0 ;-1  0 0;...
                0 0 0;1  1 0;-1 -1 0;...
                0 0 0;0 1 0 ;0 -1 0;...
                ];
            Coeffs_list = (1/(6*R_struct.b^2))*[...
                -4;2;2;...
                8;-4;-4;...
                8;-4;-4;...
                ];
        case 'x*y'
            vector_list =  [...
                1  1 0;-1 -1 0;...
                0 1 0 ;0 -1 0;...
                ];
            Coeffs_list = (1/(3*R_struct.a*R_struct.a))*[
                -sqrt(3);-sqrt(3);...
                sqrt(3);sqrt(3);...
                ];
        case 'Z'
            vector_list =  [0 0 0;0 0 1;0 0 -1];
            Coeffs_list = (1/R_struct.c^2)*[2;-1;-1];
        otherwise
            if HK.strcontain(H_str,"x*y") && ~strcmp(H_str,'x*y')
                H_str_tmp = strrep(H_str,'x*y','1');
                [vector_list,Coeffs_list] = HK.Hstr_mapping("x*y",R_struct,mode);
                [VL2,CL2] = HK.Hstr_mapping(H_str_tmp,R_struct,mode);
                [vector_list,Coeffs_list] =  HK.VLCL_ltimes(vector_list,Coeffs_list,VL2,CL2);
                return;
            else
                H_str_tmp = string(simplify(str2sym(H_str)));
                symvar_list = symvar(str2sym(H_str));
                if isempty(symvar_list)
                    vector_list = [0,0,0];
                    Coeffs_list = 1;
                elseif length(symvar_list) == 1
                    [vector_list,Coeffs_list] = HK.Hstr_mapping(H_str_tmp,R_struct,mode);
                else
                    [vector_list,Coeffs_list] = HK.Hstr_mapping(string(symvar_list(1)),R_struct,mode);
                    for i = 2:length(symvar_list)
                        [VL2,CL2] = HK.Hstr_mapping(string(symvar_list(i)),R_struct,mode);
                        [vector_list,Coeffs_list] = HK.VLCL_ltimes(vector_list,Coeffs_list,VL2,CL2);
                    end
                end
            end
    end
else
    disp('figure out the arbitrary kp2tb!!!');
end
end
