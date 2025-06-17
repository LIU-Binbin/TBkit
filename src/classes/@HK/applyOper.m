function [H_hk] = applyOper(H_hk,SymOper,options)
%APPLYOPER Apply symmetry operation to Hamiltonian
%
% Syntax:
%   H_hk = applyOper(H_hk,SymOper)
%   H_hk = applyOper(H_hk,SymOper,options)
%
% Inputs:
%   H_hk - HK object to transform
%   SymOper - Symmetry operation (Oper object or array)
%   options - Optional parameters:
%       generator - If true, generate full group (default=true)
%       Accuracy - Simplification tolerance (default=1e-3)
%       oneshot - Apply all operations at once (default=false)
%       center - Rotation center [x,y,z] (default=[0,0,0])
%
% Output:
%   H_hk - Transformed Hamiltonian
%
% Description:
%   Applies symmetry operations to the Hamiltonian, either individually or
%   by generating the full symmetry group. Handles both numeric and symbolic
%   coefficients automatically.
arguments
    H_hk HK;
    SymOper Oper = Oper();
    options.generator  logical = true;
    options.Accuracy = 1e-3;
    options.oneshot = false;
    options.center = [0,0,0];
end
if isequal(sym(zeros(size(H_hk.HcoeL))),H_hk.HcoeL)
    coe_label = false;
else
    coe_label = true;
end
if ~coe_label
    H_hk = H_hk.init();
    H_hk = H_hk.hermitize();
end
try
    BasisFunction = BasisFunc(H_hk);
    SymOper = SymOper.Ugen(BasisFunction,'Rm',H_hk.Rm,'center',options.center);
catch
    if ~isnan(SymOper.U)

    else
        error('Oper U is nan');
    end
end
if length(SymOper) == 1
    nSymOper = length(SymOper);
    i = 1;
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
            %fprintf('    ');
            H_hk = H_hk.applyOper(SymOper(i),'generator','false');
        end
        fprintf('----------   SymVarNum: %d   ----------\n',length(symvar(H_hk.HcoeL)));
    end
end
end
