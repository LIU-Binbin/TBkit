%==========================================================================
% applyOper - Apply symmetry operations to Hamiltonian in k-space
%==========================================================================
% This function applies a set of symmetry operations to a given Hamiltonian
% in k-space (H_hk) and returns the symmetrized Hamiltonian. It supports
% both generator-based and full group symmetry applications with progress
% tracking and numerical simplification options.
%
% Syntax:
%   [H_hk] = applyOper(H_hk, SymOper, options)
%
% Input Arguments:
%   H_hk    - HK object representing the Hamiltonian in k-space
%   SymOper - Oper object or array defining symmetry operations (default: Oper())
%   options - Name-value pairs:
%       generator  - Logical flag to use generator-based symmetry (default: true)
%       Accuracy   - Numerical tolerance for coefficient simplification (default: 1e-3)
%       oneshot    - Logical flag for single-step group application (default: false)
%       center    - Center coordinates for basis transformation (default: [0,0,0])
%
% Output Arguments:
%   H_hk    - Modified HK object after symmetry application
%
% Methods Called:
%   HK.init, HK.hermitize, HK.simplify, Oper.Ugen, Oper.generate_group
%
% Example:
%   H_sym = applyOper(H_initial, sym_ops, 'Accuracy', 1e-4);
%==========================================================================
function [H_hk] = applyOper(H_hk,SymOper,options)
    arguments
        H_hk HK;
        SymOper Oper = Oper();
        options.generator  logical = true;
        options.Accuracy = 1e-3;
        options.oneshot = false;
        options.center = [0,0,0];
    end
    
    % Check if coefficient matrix needs initialization
    if isequal(sym(zeros(size(H_hk.HcoeL))),H_hk.HcoeL)
        coe_label = false;
    else
        coe_label = true;
    end
    
    % Initialize and hermitize if needed
    if ~coe_label
        H_hk = H_hk.init();
        H_hk = H_hk.hermitize();
    end
    
    % Generate symmetry-adapted basis functions
    try
        BasisFunction = BasisFunc(H_hk);
        SymOper = SymOper.Ugen(BasisFunction,'Rm',H_hk.Rm,'center',options.center);
    catch
        error('Oper U is nan');
    end
    
    % Single symmetry operation case
    if length(SymOper) == 1
        nSymOper = length(SymOper);
        i = 1;
        fprintf('******** apply (%d/%d)symmetry ********\n',i,nSymOper);
        disp(SymOper(i));
        
        if options.generator
            % Generator-based symmetry application
            SymOper_tmp = SymOper(i).generate_group();
            nSymOper_tmp = length(SymOper_tmp);
            pb = CmdLineProgressBar('Applying Symmetry ...');
            H_hk_R = H_hk;
            
            for j = 1:nSymOper_tmp
                pb.print(j,nSymOper_tmp);
                [H_hk_R(j),H_hk] = applyRU(H_hk,SymOper_tmp(j));
            end
            
            pb.delete();
            H_hk = sum(H_hk_R)/nSymOper_tmp;
            H_hk = H_hk.simplify(options.Accuracy);
        else
            % Direct symmetry application
            [H_hk,H_hk_bk] = applyRU(H_hk,SymOper);
            Equationlist_r = real(H_hk.HcoeL - H_hk_bk.HcoeL) == 0;
            Equationlist_i = imag(H_hk.HcoeL - H_hk_bk.HcoeL) == 0;
            Equationlist_r = HK.isolateAll(Equationlist_r);
            Equationlist_i = HK.isolateAll(Equationlist_i);
            
            HcoeLtmp = H_hk.HcoeL ;
            HcoeLtmp_r = subs(real(HcoeLtmp),lhs(Equationlist_r),rhs(Equationlist_r));
            HcoeLtmp_i = subs(imag(HcoeLtmp),lhs(Equationlist_i),rhs(Equationlist_i));
            H_hk.HcoeL = HcoeLtmp_r + 1i*HcoeLtmp_i;
        end
        
        fprintf('----------   SymVarNum: %d   ----------\n',length(symvar(H_hk.HcoeL)));
    else
        % Multiple symmetry operations case
        if options.oneshot
            SymOper_tmp = SymOper.generate_group();
            nSymOper_tmp = length(SymOper_tmp);
            pb = CmdLineProgressBar('Applying Symmetry ...');
            H_hk_R = H_hk;
            
            for j = 1:nSymOper_tmp
                pb.print(j,nSymOper_tmp);
                [H_hk_R(j),H_hk] = applyRU(H_hk,SymOper_tmp(j));
            end
            
            pb.delete();
            H_hk = sum(H_hk_R)/nSymOper_tmp;
            H_hk = H_hk.simplify(options.Accuracy);
            return;
        end
        
        % Sequential symmetry application
        nSymOper = length(SymOper);
        for i = 1:nSymOper
            fprintf('******** apply (%d/%d)symmetry ********\n',i,nSymOper);
            disp(SymOper(i));
            
            if options.generator
                SymOper_tmp = SymOper(i).generate_group();
                nSymOper_tmp = length(SymOper_tmp);
                pb = CmdLineProgressBar('Applying Symmetry ...');
                H_hk_R = H_hk;
                
                for j = 1:nSymOper_tmp
                    pb.print(j,nSymOper_tmp);
                    [H_hk_R(j),H_hk] = applyRU(H_hk,SymOper_tmp(j));
                end
                
                pb.delete();
                H_hk = sum(H_hk_R)/nSymOper_tmp;
                H_hk = H_hk.simplify(options.Accuracy);
            else
                H_hk = H_hk.applyOper(SymOper(i),'generator','false');
            end
            
            fprintf('----------   SymVarNum: %d   ----------\n',length(symvar(H_hk.HcoeL)));
        end
    end
end

%==========================================================================
% HK Class: Hamiltonian in k-space representation
%==========================================================================
% This class represents a Hamiltonian in reciprocal space with methods for
% symmetry operations and algebraic manipulation.
%
% Properties:
%   HcoeL    - Symbolic coefficient matrix
%   Rm       - Lattice vectors matrix
%
% Methods:
%   init      - Initialize coefficient matrix
%   hermitize - Enforce Hermiticity
%   simplify  - Numerical simplification
%   applyOper - Apply symmetry operation (recursive call)
%==========================================================================

%==========================================================================
% Oper Class: Symmetry operation representation
%==========================================================================
% This class represents symmetry operations with methods for group generation
% and basis transformation.
%
% Methods:
%   Ugen          - Generate symmetry-adapted basis
%   generate_group - Generate full symmetry group from generators
%==========================================================================

%==========================================================================
% [HK Methods]
%==========================================================================

%--------------------------------------------------------------------------
% HK.init - Initialize coefficient matrix
%--------------------------------------------------------------------------
% Initializes Hamiltonian coefficient matrix to zero
%
% Syntax:
%   H_hk = init(H_hk)
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% HK.hermitize - Enforce Hermiticity
%--------------------------------------------------------------------------
% Enforces Hermitian property on the Hamiltonian
%
% Syntax:
%   H_hk = hermitize(H_hk)
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% HK.simplify - Numerically simplify coefficients
%--------------------------------------------------------------------------
% Simplifies coefficients using specified numerical tolerance
%
% Syntax:
%   H_hk = simplify(H_hk, tolerance)
%--------------------------------------------------------------------------

%==========================================================================
% [Oper Methods]
%==========================================================================

%--------------------------------------------------------------------------
% Oper.Ugen - Generate symmetry-adapted basis
%--------------------------------------------------------------------------
% Creates basis functions compatible with symmetry operations
%
% Syntax:
%   SymOper = Ugen(BasisFunction, Name, Value)
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Oper.generate_group - Generate symmetry group
%--------------------------------------------------------------------------
% Generates full symmetry group from generator elements
%
% Syntax:
%   SymGroup = generate_group()
%--------------------------------------------------------------------------
