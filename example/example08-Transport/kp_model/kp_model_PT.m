%% useful_tools
sigma_0 = pauli_matrix(0); sigma_x = pauli_matrix(1); sigma_y = pauli_matrix(2); sigma_z = pauli_matrix(3);
tau_0   = pauli_matrix(0); tau_x   = pauli_matrix(1); tau_y   = pauli_matrix(2); tau_z   = pauli_matrix(3);
%%
syms k_x k_y k_z real
syms v1 v2 v3 w delta real 
%% PT-symmetric model
H_kp = HK(4,1);
H_kp = H_kp...
     + Term( w*k_x, tau_0*sigma_0)...
     + Term(v1*k_x, tau_x*sigma_0)...
     + Term(v2*k_y, tau_y*sigma_x)...
     + Term(delta , tau_z*sigma_0);
%%
v1 = 1e6 *constants.hbar_eV_s*1e10;
v2 = 1e6 *constants.hbar_eV_s*1e10;
v3 = 0;
delta = 40e-3;
w = 0.4*v1;
%%
H_kp_n = H_kp.Subsall();
H_kp_n.bandplot([-0.2 0.2]);
%% conductivity of the second-order anomalous Hall effect
if 1 == 0
    krange = 0.02; % just around the Gamma point
    kcube_bulk = krange .* [-0.5 -0.5 -0.5; 1 0 0; 0 1 0; 0 0 1];
    Nk1 = 500;
    Nk2 = Nk1;
    klist_cart = kmeshgen(H_kp_n.Rm, kcube_bulk, "Nk1", Nk1, "Nk2", Nk2, "full_edge",true);
    
    mu_list = linspace(-0.2, 0.2, 201);
    
    chi_mu = SOAHC_int(H_kp_n, [1 2 2], klist_cart, mu_list, 'ncore',4 , 'T',20);
    kcube_ratio = krange^2;
    chi_mu = chi_mu .* kcube_ratio;
    chi_mu = chi_mu .* H_kp_n.Rm(3,3) .* 0.1; % 3D to 2D, Ang to nm
    [fig, ax] = Figs(1,1,"FigSize","wide");
    plot(mu_list, chi_mu, 'LineWidth', 2, 'Color', 'r')
    xlabel("\mu (eV)")
    ylabel("\chi^{int}_{xyy} (nmAV^{-2})")
    title(['N:',num2str(Nk1),'*',num2str(Nk2)])
end
%% distribution of Berry-connection polarizability (BCP)
if 1 == 1
    krange = 0.0637*2; % just around the Gamma point
    kcube_bulk = krange .* [-0.5 -0.5 0; 1 0 0; 0 1 0; 0 0 1];
    Nk1 = 300;
    Nk2 = Nk1;
    [klist_cart,~,~,~,~,Grid] = kmeshgen(H_kp_n.Rm, ...
        "Nk1",Nk1, ...
        "Nk2",Nk2, ...,
        "kstart",kcube_bulk(1,:), ...
        "kdir1", kcube_bulk(2,:), ...
        "kdir2", kcube_bulk(3,:), ...
        "full_edge",true);
    
    % $G_{b c, n}=2 \operatorname{Re} \sum_{\ell \neq n} \frac{\left\langle v_b\right\rangle_{n \ell}\left\langle v_c\right\rangle_{\ell n}}{\left(\varepsilon_n-\varepsilon_{\ell}\right)^3}$
    BCP_list = BCP(H_kp_n, [1 1], klist_cart, 'ncore',1 );
    %
    BCP_Grid = reshape(BCP_list(:,1), [Nk1,Nk2]);
    
    BCplot2D((BCP_Grid),Grid,double(H_kp_n.Rm),'ColorCut',0.001,'shading',false);
    % figure;
    % surf(Grid(:,:,1),Grid(:,:,2),BCP_n,BCP_n,'EdgeColor','none');
    view(0,90)
    % scatter(klist_cart(:,1),klist_cart(:,2),1,BCP_n,'filled');
    
    xlabel("k_x");
    ylabel("k_y");
    title(['G_{xx}, N:',num2str(Nk1),'*',num2str(Nk2)]);
    xlim(krange*[-0.5,0.5]);
    ylim(krange*[-0.5,0.5]);
    
    figure();
    syms k_x k_y k_z real;
    k = sqrt(k_x.^2+k_y.^2+k_z.^2);
    
    Gxx = -1/4 * (k.^2.*k_y.^2 + k_x.^2.*k_z.^2)/(k.^5.*(k_x.^2+k_y.^2));
    Gxx_0 = subs(Gxx,k_z,0);
    fsurf(Gxx_0,0.01*krange*[-0.5,0.5,-0.5,0.5],'EdgeColor','none');
    view(0,90);
    colormap("cool")

end
