function F = BerryPhaseLine_fun(Fun,kloopr,options)
%BERRYPHASELINE_FUN Calculate Berry phase along a k-loop numerically
%
%   Syntax:
%       F = BerryPhaseLine_fun(Fun,kloopr,options)
%
%   Description:
%       Computes Berry phase along a closed k-path using either:
%       1) Wavefunction overlaps (Wfun mode) or
%       2) Precomputed Berry connection (Ffun mode)
%
%   Inputs:
%       Fun      - Function handle (wavefunction or Berry connection)
%       kloopr   - Array of k-points defining the loop [N×3]
%       options  - Options structure:
%                  Rm        - Real-space lattice vectors
%                  oneshot   - Plot all at once flag (default=true)
%                  plot      - Visualization flag (default=false)
%                  car       - Cartesian coordinates flag (default=true)
%                  funtype   - Function type 'Wfun' or 'Ffun'
%
%   Output:
%       F - Berry phase (mod 2π)
arguments
    Fun function_handle;
    kloopr double;
    options.Rm = [];
    options.oneshot = true;
    options.plot = false;
    options.car = true;
    options.funtype = 'Wfun';
end
if strcmp(options.funtype,'Wfun')
    WfunType = true;
    Wfun = Fun;
    FfunType = false;
elseif strcmp(options.funtype,'Ffun')
    FfunType = true;
    Ffun = Fun;
    WfunType=false;
end
if isempty(options.Rm)
    Rm = POSCAR_read;
else
    Rm = options.Rm;
end
if options.car
else
    Gk = 2*pi*eye(3)/Rm;
    kloopr = kloopr * Gk;
end
if size(kloopr,2)<3
    kloopr = [kloopr,zeros(size(kloopr,1),3-size(kloopr,2))];
end
kn = size(kloopr,1);
if WfunType
    % $\phi \equiv-\operatorname{Im} \ln \left[\left\langle u_{0} \mid u_{1}\right\rangle\left\langle u_{1} \mid u_{2}\right\rangle \cdots\left\langle u_{N-1} \mid u_{0}\right\rangle\right]=-\sum_{j=0}^{N-1} \operatorname{Im} \ln \left\langle u_{j} \mid u_{j+1}\right\rangle$
    W1 = Wfun(kloopr(1,1),kloopr(1,2),kloopr(1,3));
    WL = zeros(size(W1,1),kn);
    for i = 1:kn
        WL(:,i) =  Wfun(kloopr(i,1),kloopr(i,2),kloopr(i,3));
    end
    dF = -imag(log(TBkit.BerryConnection(WL(:,1),WL(:,2))));
    F = dF;
elseif FfunType
    dkloopr = diff(kloopr);
    %dkloopr =[kloopr(1,:) - kloopr(end,:);dkloopr];
    dF = Ffun(kloopr(1,1),kloopr(1,2),kloopr(1,3),dkloopr(1,1),dkloopr(1,2),dkloopr(1,3));
    F = dF;
end
warning off;
if options.plot
    dBF = sum(dF);
    fig = figure('PaperType','a4letter','PaperSize',[16 8],'Color','white','Units','normalized','Position',[0.1,0.1,0.8,0.6]);
    ax0 = subplot(1,3,1,'LineWidth',1.5,'FontSize',24,'FontName',"Helvetica",'Parent',fig);
    title(ax0,'integral path');
    box(ax0,'on');
    hold(ax0,'all');
    axis(ax0,'equal');
    [fig,ax0] = BZplot(Rm,'r',0.2,fig,ax0);
    ax1 = subplot(1,3,2,'LineWidth',1.5,'FontSize',24,'FontName',"Helvetica",'Parent',fig);
    hold(ax1,'all');
    xlabel(ax1,'kpath num');
    ylabel(ax1,'dBerryphase');
    ax2 = subplot(1,3,3,'LineWidth',1.5,'FontSize',24,'FontName',"Helvetica",'Parent',fig);
    hold(ax2,'all');
    xlabel(ax2,'kpath num');
    ylabel(ax2,'Berryphase');
    count = 0;
    dBf_old = dBF;
    K_last = kloopr(1,:);
    axis(ax0,'equal');
end
for i = 2:kn-1
    k_x = kloopr(i,1);
    k_y = kloopr(i,2);
    k_z = kloopr(i,3);
    K = kloopr(i,:);
    dk = norm(K-K_last); % test
    if WfunType
        dF = -imag(log(TBkit.BerryConnection(WL(:,i),WL(:,i+1))))*dk;
    elseif FfunType
        dk_x = dkloopr(i,1);
        dk_y = dkloopr(i,2);
        dk_z = dkloopr(i,3);
        dF = Ffun(k_x,k_y,k_z,dk_x,dk_y,dk_z)*dk;
    end
    F = F +dF;
    if options.plot
        dBF = real(sum(dF));
        count = count +1;
        line(ax0,[K_last(1),k_x],[K_last(2),k_y],[K_last(3),k_z],'LineWidth',3);
        %(ax0,[K_last(1),k_x],[K_last(2),k_y],[K_last(3),k_z],'LineWidth',3);
        K_last = K;
        xlabel(ax1,string(K(1))+","+string(K(2))+","+string(K(3)));
        area(ax1,[count-1,count],[dBf_old,dBF]);
        dBf_old= dBF;
        stem(ax2,count,sum(F));
        if options.oneshot
        else
            drawnow;
        end
    end
end
if options.plot
    title(ax2,'ChernNumber = '+string((mod(real(sum(F)),(2*pi)))));
end
F = mod(real(F),2*pi);
end