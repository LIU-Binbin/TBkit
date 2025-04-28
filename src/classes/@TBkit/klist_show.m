function varargout = klist_show(TBkitobj,options)
            arguments
                TBkitobj vasplib;
                options.klist double = TBkitobj.klist_cart;
                options.ax = handle([]);
                options.show = true;
                options.color = [rand rand rand];
                options.LineWidth = 2;
            end
            if options.show
                if isempty(options.ax)
                    [~,ax] = vasplib.BZplot(TBkitobj.Rm);
                else
                    ax = options.ax;
                end
            end
            if options.show
                plot3(options.klist(:,1),options.klist(:,2),options.klist(:,3),...
                    'Color',options.color,...
                    'DisplayName','k-path',...
                    'LineWidth',options.LineWidth);
                try % test
                    for i = 1:length(TBkitobj.kpoints_name)-1
                        kpoints_r = TBkitobj.kpoints_frac(2*i-1,:)*TBkitobj.Gk;
                        text(ax,kpoints_r(1),kpoints_r(2),kpoints_r(3),TBkitobj.kpoints_name(i),'FontSize',24);
                    end
                    kpoints_r = TBkitobj.kpoints_frac(end,:)*TBkitobj.Gk;
                    text(ax,kpoints_r(1),kpoints_r(2),kpoints_r(3),TBkitobj.kpoints_name(end),'FontSize',24);
                catch

                end
                varargout{1} = ax;
            end
        end

