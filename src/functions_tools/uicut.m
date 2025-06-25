function disk_cell = uicut(disk_cell,Car_center,options)
% SIMPLEAPP Interactively explore plotting functions
%   Choose the function used to plot the sample data to see the
%   differences between surface plots, mesh plots, and waterfall plots
arguments
    disk_cell
    Car_center
    options.nside = 6;
    options.SideLength = 22.5;
    options.ShapeCenter = [];
    options.Angle = 0;
    options.ColorOrder = [1 0 0;0 0 1;0 1 0];
    options.SecondSideLength = 0;
    options.ParallelogramAngle = 120;
    options.RotateAngle = 0;
    options.RotateLayer = 1;
    options.SquareMode = 0;
    options.pgon = [];
    options.RotateMode = 1;
    options.Title = "";
    options.no_ui = false;
    options.POSCARgen = true;
    % options.twist_struct = ;
end


%%
if isa(disk_cell,'cell')
    mulitilayer_mode  = true;
else
    mulitilayer_mode  = false;
end

%
% global gldata;
% global a_handle;

%%
ColorOrder = options.ColorOrder ;
nside = options.nside;
SideLength = options.SideLength;
if options.SecondSideLength == 0
    SecondSideLength = SideLength;
else
    SecondSideLength = options.SecondSideLength ;
end
if isempty(options.ShapeCenter )
    ShapeCenter = Car_center;
else
    ShapeCenter = options.ShapeCenter ;
end
Angle = options.Angle;
ParallelogramAngle = options.ParallelogramAngle;
SquareMode = options.SquareMode;
RotateMode = options.RotateMode;
RotateAngle = options.RotateAngle;
RotateLayer = options.RotateLayer;

if options.no_ui
    orbL = rotatelayer();
    [pgon2,xv,yv] = pgon();
    if mulitilayer_mode
        for i = 1:length(disk_cell)
            Car_orbL = orbL{i}*disk_cell{i}.Rm;
            [SelectL1,SelectL2] = inpolygon(Car_orbL(:,1),Car_orbL(:,2),xv,yv);
            orbL{i} = orbL{i}(SelectL1|SelectL2,:);
        end
    else
        Car_orbL = orbL{1}*disk_cell.Rm;
        [SelectL1,SelectL2] = inpolygon(Car_orbL(:,1),Car_orbL(:,2),xv,yv);
        orbL{1} = orbL{1}(SelectL1|SelectL2,:);
        % Rm = disk_cell.Rm;
    end
    if mulitilayer_mode
        for i = 1:length(disk_cell)
            Car_orbL = orbL{i}*disk_cell{i}.Rm;
            [SelectL1,SelectL2] = inpolygon(Car_orbL(:,1),Car_orbL(:,2),xv,yv);
            disk_cell{i} = disk_cell{i}.reseq(SelectL1|SelectL2);
            disk_cell{i}.orbL = orbL{i};
        end
    else
        Car_orbL = orbL{1}*disk_cell.Rm;
        [SelectL1,SelectL2] = inpolygon(Car_orbL(:,1),Car_orbL(:,2),xv,yv);
        disk_cell = disk_cell.reseq(SelectL1|SelectL2);
        disk_cell.orbL = orbL{1};
        disk_final = disk_cell;
        % Rm = disk_cell.Rm;
    end
    if mulitilayer_mode
        disk_final = disk_cell{1};
        disk_length = length(disk_cell);
        wannum = disk_final.WAN_NUM;
        for i = 2:disk_length
            wannum = wannum + disk_cell{i}.WAN_NUM;
            disk_final.orbL = [disk_final.orbL;disk_cell{i}.orbL];
            disk_final.quantumL = [disk_final.quantumL;disk_cell{i}.quantumL]; % update orbL
            disk_final.elementL = [disk_final.elementL;disk_cell{i}.elementL]; % update orbL
        end
        disk_final.coe = false;
        disk_final.num = true;
        disk_final.HnumL = zeros(wannum);
    else
        save('disk_cell.mat','disk_cell','-v7.3');
        % fig = uifigure;
    end
    if options.POSCARgen
        disk_final.POSCAR_gen('POSCAR_bilayer.vasp');
    end
    %uialert(fig,"Success!","Alert");
    disp("Success!");
    return;
end
%%
% Create figure window
fig = uifigure('Position',[100,100,1920,1080]);
fig.CloseRequestFcn = @(src,event)my_closereq(src);
fig.Name = options.Title;
fig.WindowStyle = 'normal';
%$ Manage app layout
LayoutRow = 5;
LayoutColumn = 4;
gl = uigridlayout(fig,[LayoutRow,LayoutColumn]);
gl.RowHeight = {22,22,22,22,'1x'};
gl.ColumnWidth = {'1x','1x','1x','1x'};

% Position axes
ax = uiaxes(gl);
ax.Layout.Row = LayoutRow;
ax.Layout.Column = [1 LayoutColumn];
hold(ax,'on');

% Create UI components btn1
btn1 = uibutton(gl,"ButtonPushedFcn", @(src,event) plotButtonPushed1(src,event));
btn1.Text = 'PlotArea';
btn1.Layout.Row = 1;
btn1.Layout.Column = LayoutColumn;

% Create UI components btn2
btn2 = uibutton(gl,"ButtonPushedFcn", @(src,event) plotButtonPushed2(src,event));
btn2.Text = 'Cut';
btn2.Layout.Row = 2;
btn2.Layout.Column = LayoutColumn;

% Create UI components btn3
btn3 = uibutton(gl,"ButtonPushedFcn", @(src,event) plotButtonPushed3(src,event));
btn3.Text = 'Generate';
btn3.Layout.Row = 3;
btn3.Layout.Column = LayoutColumn;

% Create UI components btn4
btn4 = uibutton(gl,"ButtonPushedFcn", @(src,event) plotButtonPushed4(src,event));
btn4.Text = 'Rotate';
btn4.Layout.Row = 4;
btn4.Layout.Column = LayoutColumn;

% Create UI components cbx1
cbx1 = uicheckbox(gl,"ValueChangedFcn", @(src,event) checkBoxChanged1(src,event));
cbx1.Value = SquareMode;
cbx1.Text = 'SquareMode';
cbx1.Layout.Row = 3;
cbx1.Layout.Column = 1;
% Create UI components cbx2
cbx2 = uicheckbox(gl,"ValueChangedFcn", @(src,event) checkBoxChanged2(src,event));
cbx2.Value = RotateMode;
cbx2.Text = 'RotateMode';
cbx2.Layout.Row = 4;
cbx2.Layout.Column = 1;
%
% Field_infomation_struct;
maxField = 9;
maxLayer = length(disk_cell);

Field_infomation_struct = struct();
Field_infomation_struct(1).Type = "numeric";
Field_infomation_struct(1).Value = Angle;
Field_infomation_struct(1).ValueDisplayFormat = "Angle: %.2f degree";
Field_infomation_struct(1).Limits = [-360 360];
Field_infomation_struct(1).Row = 1;
Field_infomation_struct(1).Column = 1;

Field_infomation_struct = repmat(Field_infomation_struct(1),[maxField 1]);
% Create UI components Field2
Field_infomation_struct(2).Value = nside;
Field_infomation_struct(2).ValueDisplayFormat = "nside: %d ";
Field_infomation_struct(2).Limits = [3 Inf];
Field_infomation_struct(2).Row = 1;
Field_infomation_struct(2).Column = 2;

% Create UI components Field3
Field_infomation_struct(3).Value = SideLength;
Field_infomation_struct(3).ValueDisplayFormat = "SideLength: %.2f Angstrom ";
Field_infomation_struct(3).Limits = [0 Inf];
Field_infomation_struct(3).Row = 1;
Field_infomation_struct(3).Column = 3;

% Create UI components Field4
Field_infomation_struct(4).Value = ShapeCenter(1);
Field_infomation_struct(4).ValueDisplayFormat = "Center/Start(x): %.2f ";
Field_infomation_struct(4).Limits = [-Inf Inf];
Field_infomation_struct(4).Row = 2;
Field_infomation_struct(4).Column = 1;

% Create UI components Field5
Field_infomation_struct(5).Value = ShapeCenter(2);
Field_infomation_struct(5).ValueDisplayFormat = "Center/Start(y): %.2f ";
Field_infomation_struct(5).Limits = [-Inf Inf];
Field_infomation_struct(5).Row = 2;
Field_infomation_struct(5).Column = 2;

% Create UI components Field6
Field_infomation_struct(6).Value = ParallelogramAngle;
Field_infomation_struct(6).ValueDisplayFormat = "Parallelogram Angle: %.2f degree";
Field_infomation_struct(6).Limits = [-360 360];
Field_infomation_struct(6).Row = 3;
Field_infomation_struct(6).Column = 2;

% Create UI components Field7
Field_infomation_struct(7).Value = SecondSideLength;
Field_infomation_struct(7).ValueDisplayFormat = "2nd SideLength: %.2f Angstrom";
Field_infomation_struct(7).Limits = [-Inf Inf];
Field_infomation_struct(7).Row = 3;
Field_infomation_struct(7).Column = 3;

% Create UI components Field8
% Field_infomation_struct(8).Type = "text";
Field_infomation_struct(8).Value = (RotateAngle(1));
Field_infomation_struct(8).ValueDisplayFormat = "RotateAngle: %.2f degree";
Field_infomation_struct(8).Row = 4;
Field_infomation_struct(8).Column = 3;

% Create UI components Field9
Field_infomation_struct(9).Value = RotateLayer;
Field_infomation_struct(9).ValueDisplayFormat = "RotateLayer: %d";
Field_infomation_struct(9).Limits = [1 maxLayer];
Field_infomation_struct(9).Row = 4;
Field_infomation_struct(9).Column = 2;

% Create UI components Fields
for i = 1:maxField
    %bug Field_infomation_struct(i).Type
    % "Tooltip",Field_infomation_struct(i).ValueDisplayFormat, ...
    if strcmp(Field_infomation_struct(i).Type,"numeric")
        Field(i) = uieditfield(gl,Field_infomation_struct(i).Type , ...
            "Value",Field_infomation_struct(i).Value,...
            "ValueDisplayFormat",Field_infomation_struct(i).ValueDisplayFormat, ...
            "Limits",Field_infomation_struct(i).Limits, ...
            "LowerLimitInclusive","on", ...
            "UpperLimitInclusive","on", ...
            "ValueChangedFcn", @(src,event) editFieldValueChanged(src,event,i));
    else
        Field(i) = uieditfield(gl,Field_infomation_struct(i).Type , ...
            "HorizontalAlignment",'right',...
            "Value",Field_infomation_struct(i).Value,...
            "ValueChangedFcn", @(src,event) editFieldValueChanged(src,event,i));
    end
    Field(i).Layout.Row = Field_infomation_struct(i).Row;
    Field(i).Layout.Column = Field_infomation_struct(i).Column;
end
%%
orbL = rotatelayer();
plotwave();
[pgon2,xv,yv] = pgon();
a_handle = plot(ax,pgon2,'FaceAlpha',0.5);
axis(ax,'tight');
axis(ax,'equal');

    % Configure UI component appearance
    function plotButtonPushed1(src,event)
        plotpgon();
    end
    % Configure UI component appearance
    function plotButtonPushed2(src,event)
        cut();
    end
    % Configure UI component appearance
    function plotButtonPushed3(src,event)
        generate()
    end 
    % Configure UI component appearance
    function plotButtonPushed4(src,event)
        orbL = rotatelayer();
        plotwave();
        plotpgon();
    end
    % Configure UI component appearance
    function checkBoxChanged1(src,event)
        val = event.Value;
        SquareMode = val;
    end
    % Configure UI component appearance
    function checkBoxChanged2(src,event)
        val = event.Value;
        RotateMode = val;
    end
    function editFieldValueChanged(src,event,i)
        switch i
            case 1
                Angle = src.Value;
            case 2
                nside = src.Value;
            case 3
                SideLength = src.Value;
            case 4
                ShapeCenter(1) = src.Value;
            case 5
                ShapeCenter(2) = src.Value;
            case 6
                ParallelogramAngle = src.Value;
            case 7
                SecondSideLength = src.Value;
            case 8
                RotateAngle(RotateLayer) = (src.Value);
            case 9
                RotateLayer = src.Value;
            otherwise
        end
        
    end



    function my_closereq(fig)
        selection = uiconfirm(fig,'Close the figure window?',...
            'Confirmation');

        switch selection
            case 'OK'
                delete(fig)
            case 'Cancel'
                return
        end
    end
    function [pgon2,xv,yv] = pgon()
        if ~isempty(options.pgon)
            pgon2 = options.pgon;
        elseif SquareMode
            SecondSideLength_BiasX = SecondSideLength*cosd(ParallelogramAngle);
            SecondSideLength_BiasY = SecondSideLength*sind(ParallelogramAngle);
            XV = [ShapeCenter(1),...
                ShapeCenter(1)+SideLength,...
                ShapeCenter(1)+SideLength+SecondSideLength_BiasX,...
                ShapeCenter(1) + SecondSideLength_BiasX ...
                ];
            YV = [ShapeCenter(2),...           
                ShapeCenter(2),...
                ShapeCenter(2)+SecondSideLength_BiasY,...
                ShapeCenter(2)+SecondSideLength_BiasY,...
                ];
            pgon2 = polyshape(XV,YV);
        else
            pgon2 = nsidedpoly(nside,'Center',ShapeCenter([1,2]),'SideLength',SideLength);
        end
        pgon2 = rotate(pgon2,Angle,ShapeCenter([1,2]));
        xv = pgon2.Vertices(:,1);
        yv = pgon2.Vertices(:,2);
    end
    function plotwave()
        if mulitilayer_mode
            Rm = disk_cell{1}.Rm;
        else
            Rm = disk_cell.Rm;
        end
        orbL_length = length(orbL);
        cla(ax);
        hold(ax,'on');
        for tmpi = 1:orbL_length
            waveplot(orbL{tmpi},'Rm',Rm, ...
                'OrbColor',ColorOrder(tmpi,:),'OrbSize',3,...
                'ax',ax);
        end
        scatter3(ax,Car_center(1),Car_center(2),Car_center(3),30,'magenta','filled','h');
        view(ax,0,90);
        axis(ax,'tight');
        axis(ax,'equal');
    end
    function plotpgon()
        [pgon2] = pgon();
        try
            delete(a_handle);
        catch
        end
        a_handle = plot(ax,pgon2,'FaceAlpha',0.8);
    end
    function generate()
        orbL = rotatelayer();
        if mulitilayer_mode
            for tmpi = 1:length(disk_cell)
                Car_orbL = orbL{tmpi}*disk_cell{tmpi}.Rm;
                [SelectL1,SelectL2] = inpolygon(Car_orbL(:,1),Car_orbL(:,2),xv,yv);
                orbL{tmpi} = orbL{tmpi}(SelectL1|SelectL2,:);
            end
        else
            Car_orbL = orbL{1}*disk_cell.Rm;
            [SelectL1,SelectL2] = inpolygon(Car_orbL(:,1),Car_orbL(:,2),xv,yv);
            SelectL = SelectL1|SelectL2;
            orbL{1} = orbL{1}(SelectL,:);
            % Rm = disk_cell.Rm;
        end
        plotwave();
        [pgon2,xv,yv] = pgon();
        a_handle = plot(ax,pgon2,'FaceAlpha',0.5);
        axis(ax,'tight');
        axis(ax,'equal');
        
        if mulitilayer_mode
            for tmpi = 1:length(disk_cell)
                disk_cell{tmpi} = disk_cell{tmpi}.reseq(SelectL);
                disk_cell{tmpi}.orbL = orbL{tmpi};
            end
        else
            disk_cell = disk_cell.reseq(SelectL);
            disk_cell.orbL = orbL{1};
            % Rm = disk_cell.Rm;
        end
        if mulitilayer_mode
            disk_final = disk_cell{1};
            disk_length = length(disk_cell);
            wannum = disk_final.WAN_NUM;
            for tmpi = 2:disk_length
                wannum = wannum + disk_cell{tmpi}.WAN_NUM;
                disk_final.orbL = [disk_final.orbL;disk_cell{tmpi}.orbL];
                disk_final.quantumL = [disk_final.quantumL;disk_cell{tmpi}.quantumL]; % update orbL
                disk_final.elementL = [disk_final.elementL;disk_cell{tmpi}.elementL]; % update orbL
            end
            disk_final.coe = false;
            disk_final.num = true;
            disk_final.HnumL = zeros(wannum);
        else
            save('disk_cell.mat','disk_cell','-v7.3');
            disk_final = disk_cell;
            % fig = uifigure;
        end
        if options.POSCARgen
            disk_final.POSCAR_gen('POSCAR_bilayer.vasp');
        end
         uialert(fig,"Success!","Alert");
    end
    function cut()
        orbL = rotatelayer();
        if mulitilayer_mode
            for tmpi = 1:length(disk_cell)
                Car_orbL = orbL{tmpi}*disk_cell{tmpi}.Rm;
                [SelectL1,SelectL2] = inpolygon(Car_orbL(:,1),Car_orbL(:,2),xv,yv);
                orbL{tmpi} = orbL{tmpi}(SelectL1|SelectL2,:);
            end
        else
            Car_orbL = orbL{1}*disk_cell.Rm;
            [SelectL1,SelectL2] = inpolygon(Car_orbL(:,1),Car_orbL(:,2),xv,yv);
            orbL{1} = orbL{1}(SelectL1|SelectL2,:);
            % Rm = disk_cell.Rm;
        end
        plotwave();
        [pgon2,xv,yv] = pgon();
        a_handle = plot(ax,pgon2,'FaceAlpha',0.5);
        axis(ax,'tight');
        axis(ax,'equal');
    end
    function orbL = init_orbL()
        if mulitilayer_mode
            for tmpi = 1:length(disk_cell)
                orbL{tmpi} = disk_cell{tmpi}.orbL;
            end
            if length(RotateAngle) ~= length(disk_cell)
                RotateAngle = repmat(RotateAngle,[1 length(disk_cell)]);
            end
        else
            orbL{1} = disk_cell.orbL;
            % Rm = disk_cell.Rm;
        end
    end
    function orbL = rotatelayer()
        orbL = init_orbL();
        if RotateMode
            if mulitilayer_mode
                Rm = disk_cell{1}.Rm;
                for tmpi = 1:length(disk_cell)
                    % disp(RotateAngle(tmpi));
                    orbL{tmpi} = RotateOneLayer(orbL{tmpi},Rm,Car_center,RotateAngle(tmpi));
                end
            else
                orbL{1} = RotateOneLayer(orbL{1},disk_cell.Rm,Car_center,RotateAngle);
            end
        end
    end
    function OneOrbL = RotateOneLayer(OneOrbL,Rm,Car_center,theta)
        % Car
        Car_orbL = OneOrbL * Rm - Car_center;
        % disp(theta);
        % cylindrical coordinate system
        %figure();scatter3(Car_orbL(:,1),Car_orbL(:,2),Car_orbL(:,3));view(0,90);
        [thetaL,rhoL,zL] = cart2pol(Car_orbL(:,1),Car_orbL(:,2),Car_orbL(:,3));
        thetaL = thetaL + theta/180*pi ;
        [xL,yL,zL] = pol2cart(thetaL,rhoL,zL);
        %figure();scatter3(xL,yL,zL);;view(0,90);
        CarOneLayer = [xL,yL,zL] + Car_center;
        OneOrbL = CarOneLayer/Rm;
    end
end