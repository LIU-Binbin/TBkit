function [klist_cart,klist_frac] = kloop1D(kpoint_frac,Orientation,radius,opt)
            arguments
                kpoint_frac;
                Orientation = [0,0,1];
                radius = 0.01;
                opt.nk = 100;
                opt.Gk = [];
                opt.inputCar logical= false
                opt.enforceCar logical= false
                opt.theta_0 = 0;
                opt.anticlockwise = true;
            end
            if isempty(opt.Gk)
                Rm = vasplib.POSCAR_read;
                opt.Gk = (2*pi*eye(3)/Rm).';
            end
            if opt.inputCar
                kpoint_cart = kpoint_frac;
                Orientation_cart = Orientation;
                radius_cart = radius;
            else
                kpoint_cart = kpoint_frac*opt.Gk;
                Orientation_cart = Orientation*opt.Gk;
                radius_cart = norm(opt.Gk)*radius;
            end
            if opt.enforceCar
                %Orientation_r = Orientation;
            end
            %
            % centered on (x0, y0, z0), with a radius of r, and with
            % normal vector n
            n=Orientation_cart; %法向量n
            r=radius_cart; %圆的半径为1
            c=kpoint_cart; %圆心的坐标
            if opt.anticlockwise
                theta=linspace(opt.theta_0,opt.theta_0+2*pi,opt.nk)'; %theta角从0到2*pi
            else
                theta=linspace(opt.theta_0,opt.theta_0-2*pi,opt.nk)'; %theta角从0到2*pi
            end
            a=cross([0 1 0],n); %n与j叉乘，求取a向量
            if ~any(a) %如果a为零向量，将n与j叉乘
                a=cross(n,[1 0 0]);
            end
            b=cross(n,a); %求取b向量
            a=a/norm(a); %单位化a向量
            b=b/norm(b); %单位化b向量
            c1=c(1)*ones(size(theta,1),1);
            c2=c(2)*ones(size(theta,1),1);
            c3=c(3)*ones(size(theta,1),1);
            k_xL=c1+r*a(1)*cos(theta)+r*b(1)*sin(theta);%圆上各点的x坐标
            k_yL=c2+r*a(2)*cos(theta)+r*b(2)*sin(theta);%圆上各点的y坐标
            k_zL=c3+r*a(3)*cos(theta)+r*b(3)*sin(theta);%圆上各点的z坐标
            %h = plotg(k_xL,k_yL,k_zL);
            %set(h,'LineWidth',3);
            klist_cart = [k_xL,k_yL,k_zL];
            klist_frac = klist_cart/opt.Gk;
        end
