function [figure_dos, Pdensity, Totdensity] = dosplot(Ecut, choose_list, mode, options)
    % DOSPLOT Generates Density of States (DOS) and partial DOS (PDOS) plots.
    %   Inputs:
    %       Ecut: Energy range for the plot (in eV) [y-axis limits].
    %       choose_list: List of selected atoms or orbitals for PDOS.
    %       mode: Mode of operation. Default is 'vaspkit-dir'.
    %       options: Struct containing options for the plot, e.g., options.d (to include d-orbitals).
    %   Outputs:
    %       figure_dos: Handles to the created figures.
    %       Pdensity: Partial density of states data.
    %       Totdensity: Total density of states data.

    % Default argument values
    arguments
        Ecut = [-2, 2];              % Energy range for the plot (default [-2, 2] eV)
        choose_list = -1;            % Default behavior: plot all (if -1)
        mode = 'vaspkit-dir';        % Default mode: 'vaspkit-dir'
        options.d = true;            % Option to include d-orbitals (default true)
    end
    
    % Initialize necessary imports
    import park.*
    import TBkit_tool.*
    
    % Check the mode and process accordingly
    if strcmp(mode, 'vaspkit-dir')
        % Generate total DOS (TDOS)
        filename = 'tdos.dat';
        [Totdensity, EIGENCAR, ~] = WEIGHTCAR_gen(filename, -1, 'vaspkit-DOS-silence');
        
        % Generate partial DOS (PDOS) from the directory
        [Pdensity, Pdos_namelist, width] = WEIGHTCAR_gen_dir('DOS');
        
        % Determine orbital mode based on the width of the PDOS data
        if width > 10
            f_mode = true;
            num2orbital_name = orbital_maprule_gen(1);
        else
            f_mode = false;
            num2orbital_name = orbital_maprule_gen(0);
        end
        
        % If choose_list < 0, generate and plot total, element-wise, and orbital-wise DOS
        if choose_list < 0
            fprintf('Generating vaspkit DOS for the directory...\n');
            
            % Plot total DOS
            ax(1) = dos_plot_set(EIGENCAR, Totdensity, 'DOS-TOT', Ecut, 'DOS-TOT');
            figure_dos(1) = ax(1).Parent;
            
            % Prepare element-wise DOS and plot
            EIGENCAR_element(:, 1) = EIGENCAR;
            WEIGHTCAR_element(:, 1) = Totdensity;
            Name_list(1, :) = "TOT";
            
            for i = 1:length(Pdos_namelist)
                EIGENCAR_element(:, i + 1) = EIGENCAR;
                WEIGHTCAR_element(:, i + 1) = Pdensity(:, end, i);
                Name_list(i + 1, :) = Pdos_namelist(i, :);
            end
            
            [ax(2)] = dos_plot_set(EIGENCAR_element, WEIGHTCAR_element, Name_list, Ecut, 'DOS-element');
            figure_dos(2) = ax(2).Parent;
            
            % Prepare orbital-wise DOS and plot
            Name_list = [""];
            EIGENCAR_orbital(:, 1) = EIGENCAR;
            WEIGHTCAR_orbital(:, 1) = Totdensity;
            Name_list(1, :) = 'TOT';
            
            % Sum s-orbital contribution
            EIGENCAR_orbital(:, 2) = EIGENCAR;
            temp_WEIGHTCAR = [];
            for i = 1:length(Pdos_namelist)
                temp_WEIGHTCAR = [temp_WEIGHTCAR, Pdensity(:, 1, i)];
            end
            WEIGHTCAR_orbital(:, 2) = sum(temp_WEIGHTCAR, 2);
            Name_list(2, :) = 's';
            
            % Sum p-orbital contribution
            EIGENCAR_orbital(:, 3) = EIGENCAR;
            temp_WEIGHTCAR = [];
            for i = 1:length(Pdos_namelist)
                temp_WEIGHTCAR = [temp_WEIGHTCAR, Pdensity(:, 2:4, i)];
            end
            WEIGHTCAR_orbital(:, 3) = sum(temp_WEIGHTCAR, 2);
            Name_list(3, :) = 'p';
            
            % Sum d-orbital contribution (if enabled in options)
            if options.d
                EIGENCAR_orbital(:, 4) = EIGENCAR;
                temp_WEIGHTCAR = [];
                for i = 1:length(Pdos_namelist)
                    temp_WEIGHTCAR = [temp_WEIGHTCAR, Pdensity(:, 5:9, i)];
                end
                WEIGHTCAR_orbital(:, 4) = sum(temp_WEIGHTCAR, 2);
                Name_list(4, :) = 'd';
            end
            
            % Sum f-orbital contribution (if applicable)
            if width > 10
                EIGENCAR_orbital(:, 5) = EIGENCAR;
                temp_WEIGHTCAR = [];
                for i = 1:length(Pdos_namelist)
                    temp_WEIGHTCAR = [temp_WEIGHTCAR, Pdensity(:, 10:16, i)];
                end
                WEIGHTCAR_orbital(:, 5) = sum(temp_WEIGHTCAR, 2);
                Name_list(5, :) = 'f';
            end
            
            [ax(3)] = dos_plot_set(EIGENCAR_orbital, WEIGHTCAR_orbital, Name_list, Ecut, 'DOS-orbital');
            figure_dos(3) = ax(3).Parent;
            
            % Generate individual orbital plots for each atom in the dataset
            count = 4;
            if options.d
                width = 9; % for d-orbitals
            else
                width = 4; % without d-orbitals
            end
            
            for i = 1:length(Pdos_namelist)
                Name_list = [""];
                EIGENCAR_atom(:, 1) = EIGENCAR;
                WEIGHTCAR_atom(:, 1) = Totdensity;
                Name_list(1, :) = "TOT";
                
                % Add orbital contributions for each atom
                WEIGHTCAR_atom(:, 2:width + 1) = Pdensity(:, 1:width, i);
                for j = 1:width
                    EIGENCAR_atom(:, j + 1) = EIGENCAR;
                    tempstring = strcat(strcat(Pdos_namelist(i, :), '-'), num2orbital_name(j));
                    Name_list(j + 1, :) = tempstring;
                end
                
                % Plot the orbital DOS for each atom
                figure_name = strcat('DOS-orbital-', Pdos_namelist(i, :));
                [ax(count)] = dos_plot_set(EIGENCAR_atom, WEIGHTCAR_atom, Name_list, Ecut, figure_name);
                figure_dos(count) = ax(count).Parent;
                count = count + 1;
            end
            
        else
            % Handle custom user selection mode
            print_PDOS_information_in_dir(Pdos_namelist, width);
            if isa(choose_list, 'double')
                fprintf("Please select manually.\n");
            else
                % Prepare custom PDOS based on user selection
                EIGENCAR_atom(:, 1) = EIGENCAR;
                WEIGHTCAR_atom(:, 1) = Totdensity;
                Name_list(1, :) = "TOT";
                Nchoose = length(choose_list);
                
                for i = 1:Nchoose
                    EIGENCAR_atom(:, i + 1) = EIGENCAR;
                    Name_list(i + 1, :) = choose_list(i).displayname;
                    WEIGHTCAR_atom(:, i + 1) = sum(sum(Pdensity(:, ...
                        choose_list(i).obital_seq_list, ...
                        choose_list(i).atom_seq_list), 3), 2);
                end
                
                % Generate PDOS plot for the selected atoms/orbitals
                figure_name = 'PDOS';
                [ax] = dos_plot_set(EIGENCAR_atom, WEIGHTCAR_atom, Name_list, Ecut, figure_name);
                figure_dos = ax.Parent;
            end
        end
    end
end



% function [figure_k,axis] = dos_plot_one(EIGENCAR,WEIGHTCAR,Name,axis,figure_k)
%     Fontname="Helvetica";
%     if nargin < 4
%         figure_k= figure('PaperType','a4letter','PaperSize',[8 8],'Color','white');
%         axis = axes('Parent',figure_k,'LineWidth',1.5,'FontSize',24.0,'FontName',Fontname);         %Position [left bottom width height]
%     end
%     plot(axis,EIGENCAR,WEIGHTCAR,'LineWidth',2.0,'color',[rand,rand,rand],'DisplayName',Name);
% end