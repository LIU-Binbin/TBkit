function [H_hk] = applyOper(H_hk,SymOper,options)
%APPLYOPER Apply symmetry operations to Hamiltonian in k-space
%   [H_hk] = applyOper(H_hk, SymOper, options) applies symmetry operations to
%   a k-space Hamiltonian object.
%
%   Inputs:
%       H_hk    - HK object representing the Hamiltonian
%       SymOper - Symmetry operation (Oper object or array)
%       options - Optional parameters:
%           .generator  : Use group generators (default: true)
%           .Accuracy   : Simplification tolerance (default: 1e-3)
%           .oneshot    : Apply all symmetries at once (default: false)
%           .center     : Center for symmetry operations (default: [0,0,0])
%
%   Output:
%       H_hk - Symmetry-adapted Hamiltonian
%
%   Example:
%       H_sym = applyOper(H, [sym1,sym2], 'Accuracy', 1e-5);
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
