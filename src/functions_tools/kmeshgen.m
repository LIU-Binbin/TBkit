function [klist_cart, klist_frac, klist_k1, klist_k2, klist_k3, Grid] = kmeshgen(Rm, KCUBE_BULK, options)
arguments
    Rm (3,3) double = []
    KCUBE_BULK double = []    % [kstart; kdir1; kdir2; kdir3]
    options.kstart(1,3) double = [0 0 0]
    options.kdir1 (1,3) double = [1 0 0]
    options.kdir2 (1,3) double = [0 1 0]
    options.kdir3 (1,3) double = [0 0 1]
    options.Nk1 double = 1
    options.Nk2 double = 1
    options.Nk3 double = 1
    options.full_edge logical = false
    options.dimension int8 = 3
end
%%
if isempty(Rm)
    try
        Rm = POSCAR_read;
    catch
        warning("POSCAR read error, we will use an identity matrix as Rm")
        Rm = eye(3);
    end
end

if ~isempty(KCUBE_BULK)
    kstart = KCUBE_BULK(1,:);
    vk = KCUBE_BULK(2:4,:);
else
    kstart = options.kstart;
    vk = [options.kdir1; options.kdir2; options.kdir3];
end
Nkvec = [options.Nk1, options.Nk2, options.Nk3];
%%
if options.dimension == 2
    vk(3,:) = 0;
    Nkvec(3) = 1;
elseif options.dimension == 3

end
%%
if options.full_edge
    klist_k1 = linspace3([0 0 0], vk(1,:), Nkvec(1));
    klist_k2 = linspace3([0 0 0], vk(2,:), Nkvec(2));
    klist_k3 = linspace3([0 0 0], vk(3,:), Nkvec(3));
else
    klist_k1 = linspace3([0 0 0], vk(1,:), Nkvec(1)+1);
    klist_k2 = linspace3([0 0 0], vk(2,:), Nkvec(2)+1);
    klist_k3 = linspace3([0 0 0], vk(3,:), Nkvec(3)+1);
end
%%
klist_frac = zeros(Nkvec(1) * Nkvec(2) * Nkvec(3), 3);
for c = 1:Nkvec(3)
    start1 = Nkvec(1)*Nkvec(2) * (c-1);
    klist_frac( start1+1: Nkvec(1)*Nkvec(2)+start1, 1:3)=...
        klist_frac( start1+1: Nkvec(1)*Nkvec(2)+start1, 1:3) + klist_k3(c,1:3);
    for b = 1:Nkvec(2)
        start2 = start1 + (b-1)*Nkvec(1);
        klist_frac( start2+1: Nkvec(1)+start2, 1:3)=...
            klist_frac( start2+1: Nkvec(1)+start2, 1:3) +...
            klist_k1(1:Nkvec(1), 1:3) + klist_k2(b, 1:3);
    end
end
klist_frac = klist_frac + kstart;
Gk = (eye(length(Rm))*2*pi/(Rm)).';
klist_cart = klist_frac * Gk;
%% 2D Grid
if options.Nk3 == 1 || all(options.kdir3 == 0)
    Grid(:,:,1)  = reshape(klist_cart(:,1), [options.Nk1, options.Nk2]);
    Grid(:,:,2)  = reshape(klist_cart(:,2), [options.Nk1, options.Nk2]);
    Grid(:,:,3)  = reshape(klist_cart(:,3), [options.Nk1, options.Nk2]);
end