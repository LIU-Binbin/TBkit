function H_hk = Subsall(H_hk, options)
arguments
    H_hk HK
    options.mode {mustBeMember(options.mode,{'num','file','gen','sym'})}= 'num';
    options.get_dH_dk_fun logical = true
end

if exist('para.mat','file') && strcmp(options.mode,'file')
    load('para.mat');
else
end
switch options.mode
    case {'num','file','gen'}
        HcoeL_temp = subs(H_hk.HcoeL);
        symname = symvar(HcoeL_temp);
        if ~isempty(symname)
            for i = 1:length(symname)
                warning('this varible: %s is not defind',string(symname(i)));
            end
            disp(symname);
            error('please introduce the varible value above');
        else
            H_hk.HnumL = double(HcoeL_temp);
        end
        H_hk.Hk_num  = subs(H_hk.Hk_sym);
        H_hk.Hsym  = H_hk.Hk_num ;
        H_hk.num = 1;
    case 'sym'
        H_hk.HcoeL = subs(H_hk.HcoeL);
        H_hk.Trig_to_save = subs(H_hk.Trig_to_save);
end
%%
% Convert the symbolic derivatives into MATLAB function handles
syms k_x k_y k_z real
H_hk.Hfun = matlabFunction(H_hk.sym, 'Vars', [k_x, k_y, k_z]);
%%
if options.get_dH_dk_fun
    symL = H_hk.HsymL;
    
    dsym_dkx = diff(symL, k_x);  % Derivative with respect to k_x
    dsym_dky = diff(symL, k_y);  % Derivative with respect to k_y
    dsym_dkz = diff(symL, k_z);  % Derivative with respect to k_z
    
    nbands = H_hk.Nbands;
    HnumL = reshape(H_hk.HnumL, nbands^2, []);
    
    dH_dk_sym = zeros(nbands, nbands, 3);
    dH_dk_sym(:,:,1) = reshape(HnumL * dsym_dkx.', nbands, nbands);
    dH_dk_sym(:,:,2) = reshape(HnumL * dsym_dky.', nbands, nbands);
    dH_dk_sym(:,:,3) = reshape(HnumL * dsym_dkz.', nbands, nbands);
    
    % Convert the symbolic derivatives into MATLAB function handles
    H_hk.dH_dk_fun = matlabFunction(dH_dk_sym, 'Vars', [k_x, k_y, k_z]);
end
end
