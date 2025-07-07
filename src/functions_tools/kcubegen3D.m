function [klist_cart, klist_frac] = kcubegen3D(options)
arguments
    options.Rm (3,3) double = [1 0 0; 0 1 0; 0 0 1];
    options.KCUBE_BULK double = [-0.5 -0.5 -0.5; 1 0 0; 0 1 0; 0 0 1]; % [original_point; vk1; vk2; vk3] 
    options.nk double = 100
end
%%
if isscalar(options.nk)
    nk = [options.nk, options.nk, options.nk];
elseif length(options.nk) == 3
    nk = options.nk;
else
    error("nk is wrongly set! It must be a scalar or a vector like [nk1 nk2 nk3].")
end
%%
original_point = options.KCUBE_BULK(1,:);
vk = options.KCUBE_BULK(2:4,:);
klist_k1 = linspace3([0 0 0], vk(1,:), nk(1)+1);
klist_k2 = linspace3([0 0 0], vk(2,:), nk(2)+1);
klist_k3 = linspace3([0 0 0], vk(3,:), nk(3)+1);
%%
klist_frac = zeros(nk(1) * nk(2) * nk(3), 3);
for c = 1:nk(3)
    start1 = nk(1)*nk(2) * (c-1);
    klist_frac( start1+1: nk(1)*nk(2)+start1, 1:3)=...
        klist_frac( start1+1: nk(1)*nk(2)+start1, 1:3) + klist_k3(c,1:3);
    for b = 1:nk(2)
        start2 = start1 + (b-1)*nk(1);
        klist_frac( start2+1: nk(1)+start2, 1:3)=...
            klist_frac( start2+1: nk(1)+start2, 1:3) +...
            klist_k1(1:nk(1), 1:3) + klist_k2(b, 1:3);
    end
end
klist_frac = klist_frac + original_point;
Gk = (eye(length(options.Rm))*2*pi/(options.Rm)).';
klist_cart = klist_frac * Gk;
end

function vlist = linspace3(v1, v2, n)
dim = length(v1);
vlist = zeros(n, dim);
for i = 1:dim
    vlist(1:n,i) = linspace(v1(i),v2(i),n);
end
end
