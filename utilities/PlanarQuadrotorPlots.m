classdef PlanarQuadrotorPlots
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        labelFontSize (1,1) {mustBeNumeric}
        tickFontSize (1,1) {mustBeNumeric,mustBeInteger}
        legendFontSize (1,1) {mustBeNumeric,mustBeInteger}
        width (1,1) {mustBeNumeric,mustBeNonnegative}
        margin (1,1) {mustBeNumeric,mustBeNonnegative}
        alphaArea (1,1) {mustBeNumeric,mustBeNonnegative}
        lineWidth (1,1) {mustBeNumeric,mustBeNonnegative}
        lineStyle (:,:) string {mustBeMember(lineStyle,{'-','--',':','-.'})} = '-'
        xLabels (:,:) cell {mustBeText,mustBeVector(xLabels,"allow-all-empties")}
        uLabels (:,:) cell {mustBeText,mustBeVector(uLabels,"allow-all-empties")}
        mapping (:,:) cell {mustBeVector(mapping,"allow-all-empties")}
    end
    properties (Access=protected,Hidden)
        uBound = @(x,y) cellfun(@(a,b)a+b,x,y,UniformOutput=false);
        lBound = @(x,y) cellfun(@(a,b)a-b,x,y,UniformOutput=false);
    end

    methods
        function obj = PlanarQuadrotorPlots(options)
            arguments
                options.labelFontSize (1,1) {mustBeNumeric} = 9;
                options.tickFontSize (1,1) {mustBeNumeric,mustBeInteger} = 6;
                options.legendFontSize (1,1) {mustBeNumeric,mustBeInteger} = 10;
                options.lineWidth (1,1) {mustBeNumeric,mustBeNonnegative} = 1;
                options.alphaArea (1,1) {mustBeNumeric,mustBeNonnegative} = 0.3;
                options.width (1,1) {mustBeNumeric,mustBeNonnegative} = 25;
                options.margin (1,1) {mustBeNumeric,mustBeNonnegative} = 1;
                options.xLabels (:,:) cell {mustBeText,mustBeVector} = {''}
                options.uLabels (:,:) cell {mustBeText,mustBeVector} = {''}
                options.mapping (:,:) cell {mustBeVector} = {[1,2],[3,6],[4,5]}
            end
            
            % Define labels
            obj.xLabels = options.xLabels;
            obj.uLabels = options.uLabels;
            if all(cellfun(@isempty, obj.xLabels))
                obj.xLabels{1} = 'x';
                obj.xLabels{2} = 'y';
                obj.xLabels{3} = '\phi';
                obj.xLabels{4} = '{\dot x}';
                obj.xLabels{5} = '{\dot y}';
                obj.xLabels{6} = '{\dot \phi}';
            end
            if all(cellfun(@isempty, obj.uLabels))
                obj.uLabels{1} = 'u_1';
                obj.uLabels{2} = 'u_2';
            end

            % Check state to plot mapping
            obj.mapping = options.mapping;
            if max(cellfun(@(x) max(x,[],'all'), obj.mapping)) > numel(obj.xLabels)
                error('Mapping must contain valid state indices.')
            end
            
            % Define appearance
            obj.labelFontSize = options.labelFontSize;
            obj.tickFontSize = options.tickFontSize;
            obj.legendFontSize = options.legendFontSize;
            obj.lineWidth = options.lineWidth;
            obj.alphaArea = options.alphaArea;
            obj.width = options.width;
            obj.margin = options.margin;
        end

        function trajectories(obj,t,x,u,xTubes,uTubes,label,options)
            arguments (Input)
                obj (1,1) PlanarQuadrotorPlots
            end
            arguments (Input,Repeating)
                t (:,:) {mustBeNumeric,mustBeNonnegative,mustBeVector}
                x (:,:) {mustBeNumeric}
                u (:,:) {mustBeNumeric}
                xTubes (:,:) {mustBeNumeric}
                uTubes (:,:) {mustBeNumeric}
                label (1,1) string {mustBeText}
            end
            arguments (Input)
                options.sharedAxes (1,1) logical = true
                options.xLimits (:,:) cell {mustBeVector} = cell(6,1)
                options.uLimits (:,:) cell {mustBeVector} = cell(2,1)
                options.xColors (:,:) cell {mustBeVector} = cell(1,1)
                options.uColors (:,:) cell {mustBeVector} = cell(1,1)
                options.mapping (:,:) cell {mustBeVector} = obj.mapping
                options.xLabels (:,:) cell {mustBeText,mustBeVector} = obj.xLabels
                options.uLabels (:,:) cell {mustBeText,mustBeVector} = obj.uLabels
                options.xAxisLabel (1,1) string = "Time [s]"
                options.labelFontSize (1,1) {mustBeNumeric,mustBeInteger} = obj.labelFontSize;
                options.tickFontSize (1,1) {mustBeNumeric,mustBeInteger} = obj.tickFontSize;
                options.legendFontSize (1,1) {mustBeNumeric,mustBeInteger} = obj.legendFontSize;
                options.legendLineLength (1,1) {mustBeInteger,mustBeInRange(options.legendLineLength,1,30)} = 15;
                options.lineWidth (1,1) {mustBeNumeric,mustBeNonnegative} = obj.lineWidth;
                options.lineStyle (:,:) string {mustBeMember(options.lineStyle,{'-','--',':','-.'})} = ["-","--",":","-."]
                options.alphaArea (1,1) {mustBeNumeric,mustBeNonnegative} = obj.alphaArea;
                options.width (1,1) {mustBeNumeric,mustBeNonnegative} = obj.width;
                options.height (1,1) {mustBeNumeric,mustBeNonnegative} = 0.32*obj.width;
            end

            % Get number of repeating arguments
            N = length(x);

            % Define colormaps (repeating)
            if any(cellfun(@isempty, options.xColors))
                options.xColors = obj.get_colormap(@cool,N,2);
            end
            if any(cellfun(@isempty, options.uColors))
                options.uColors = obj.get_colormap(@cool,N,2);
            end

            % Define lineStyles
            options.lineStyle = num2cell(options.lineStyle);

            % Compute tubes (repeating)
            xUpper = obj.uBound(x,xTubes);
            xLower = obj.lBound(x,xTubes);
            uUpper = obj.uBound(u,uTubes);
            uLower = obj.lBound(u,uTubes);
            
            % Compute shared axis limits
            if all(cellfun(@isempty, options.xLimits))
                sharedLimits = 1.1*[min([x{:}],[],"all"), ...
                                    max([x{:}],[],"all")];
            else
                sharedLimits = 1.1*[min([options.xLimits{:}],[],"all"), ...
                                    max([options.xLimits{:}],[],"all")];
            end

            % Check state to plot mapping
            idx = options.mapping;
            if max(cellfun(@(x) max(x,[],'all'), idx)) > numel(obj.xLabels)
                error('Mapping must contain valid state indices.')
            end

            % Create tiled figure
            figure()
            tiledlayout('horizontal',TileSpacing='Tight',Padding='Compact')
            
            % Plot states
            for i = 1:numel(idx)
                nexttile; hold on; grid on;
                
                % Extract cell vectors (repeating)
                x_i = cellfun(@(X) X(idx{i},:), x, UniformOutput=false);
                xLower_i = cellfun(@(X) X(idx{i},:), xLower, UniformOutput=false);
                xUpper_i = cellfun(@(X) X(idx{i},:), xUpper, UniformOutput=false);

                % Plot states (repeating)
                obj.plot_states(t,x_i, ...
                                xLower_i, ...
                                xUpper_i, ...
                                options.xLimits(idx{i},:), ...
                                options.xColors, ...
                                options.xLabels(idx{i}), ...
                                options.lineWidth, ...
                                options.lineStyle, ...
                                options.alphaArea)

                % Define labels and axis limits
                ylabel('States', Interpreter='latex');
                xlabel(options.xAxisLabel, Interpreter='latex');
                ylim('padded');
                xlim('tight');

                % Define axes
                if options.sharedAxes
                    ylim(sharedLimits);
                end
                if i~=1 && options.sharedAxes
                    set(gca, YTickLabel=[], YColor='none');
                    set(gca().YLabel, String='')
                end

                % Define fontsizes
                set(gca().XAxis, FontSize=options.tickFontSize)
                set(gca().YAxis, FontSize=options.tickFontSize)
                set(gca().XLabel, FontSize=options.labelFontSize)
                set(gca().YLabel, FontSize=options.labelFontSize)
                
                % Plot legend
                lgd = legend(Interpreter='latex', ...
                             Location='northoutside', ...
                             Orientation='horizontal', ...
                             FontSize=options.legendFontSize, ...
                             Box='off');
                lgd.ItemTokenSize = [options.legendLineLength;0];
            end
            
            % Plot inputs (repeating)
            nexttile; hold on; grid on;
            obj.plot_inputs(t,u, ...
                            uLower, ...
                            uUpper, ...
                            options.uLimits, ...
                            options.uColors, ...
                            options.uLabels, ...
                            options.lineWidth, ...
                            options.lineStyle, ...
                            options.alphaArea)

            % Define labels and axis limits
            ylabel('Inputs', Interpreter='latex');
            xlabel(options.xAxisLabel, Interpreter='latex');
            ylim('padded');
            xlim('tight');
            
            % Define fontsizes
            set(gca().XAxis, FontSize=options.tickFontSize)
            set(gca().YAxis, FontSize=options.tickFontSize)
            set(gca().XLabel, FontSize=options.labelFontSize)
            set(gca().YLabel, FontSize=options.labelFontSize)
        
            % Plot legend
            lgd = legend(Interpreter='latex', ...
                         Location='northoutside', ...
                         Orientation='horizontal', ...
                         FontSize=options.legendFontSize, ...
                         Box='off');
            lgd.ItemTokenSize = [options.legendLineLength;0];

            % Plot invisible dummy axis on the right for legend
            nexttile; hold on;
            set(gca(), Visible='off')
            for i = 1:N
                patch(0,0,options.xColors{i}{1}, ...
                          FaceAlpha=options.alphaArea, ...
                          EdgeColor=options.xColors{i}{2}, ...
                          LineWidth=1.5*options.lineWidth, ...
                          DisplayName=label{i});
            end
            
            % Plot legend
            lgd = legend(Interpreter='latex', ...
                         Location='west', ...
                         Orientation='vertical', ...
                         FontSize=options.legendFontSize, ...
                         Box='off');
            lgd.Layout.Tile = numel(idx)+2;

            % Set plot width and plot height
            set(gcf, Units='centimeters');
            set(gcf, Position=[obj.margin obj.margin options.width options.height]);
        end

        function errors(obj,noiseErr,paramErr,linErr,options)
            arguments
                obj (1,1) PlanarQuadrotorPlots
                noiseErr (:,:) {mustBeNumeric}
                paramErr (:,:) {mustBeNumeric}
                linErr (:,:) {mustBeNumeric}
                options.colors (:,:) cell {mustBeVector} = {"#0C0786";"#9B179E";"#EC7853";"#EFF821"};
                options.labels (:,:) cell {mustBeText,mustBeVector} = obj.xLabels
                options.labelFontSize (1,1) {mustBeNumeric} = obj.labelFontSize;
                options.tickFontSize (1,1) {mustBeNumeric,mustBeInteger} = obj.tickFontSize;
                options.legendFontSize (1,1) {mustBeNumeric,mustBeInteger} = obj.legendFontSize;
                options.lineWidth (1,1) {mustBeNumeric,mustBeNonnegative} = obj.lineWidth;
                options.alphaArea (1,1) {mustBeNumeric,mustBeNonnegative} = obj.alphaArea;
                options.width (1,1) {mustBeNumeric,mustBeNonnegative} = obj.width;
                options.height (1,1) {mustBeNumeric,mustBeNonnegative} = 0.2*obj.width;
            end

            % Create tiled figure
            figure()
            tiledlayout('horizontal', TileSpacing='compact')

            % Stack errors for bar plot
            err_stacked = cat(3,noiseErr,paramErr,linErr);

            % Compute shared axis limits
            max_rounded = 10^ceil(log10(max(sum(err_stacked,3),[],'all')));
            sharedLimits = [0.0005, max_rounded];

            % Define color labels
            colorLabels = {'$\mathcal{W}$', '${\Theta}^\star$', '$\tau^{\star 2} \mu$'};
            
            % Plot stacked bars
            for i = 1:size(err_stacked,1)
                nexttile; hold on; grid on;
                b = bar(squeeze(err_stacked(i,:,:)), 'stacked');
                
                % Set colors and labels
                for j = 1:size(err_stacked,3)
                    b(j).FaceColor = options.colors{j};
                    b(j).DisplayName = colorLabels{j};
                end

                % Define labels and axis limits
                ylim(sharedLimits);
                xlabel('Step $k$', Interpreter='latex', FontSize=options.labelFontSize);
                ylabel('$\sigma_{i,k}^\star$', Interpreter='latex', FontSize=options.labelFontSize);
                title(['$',options.labels{i},'$'], Interpreter='latex', FontSize=options.labelFontSize);
                set(gca, YScale='log', XTickLabel=[], FontSize=options.tickFontSize)

                % Set yticks
                if i ~= 1
                    set(gca, YTickLabel=[], YColor='none');
                    set(gca().YLabel, String='')
                end
            end
            
            % Plot legend
            lgd = legend(Interpreter='latex', ...
                         FontSize=options.legendFontSize, ...
                         Box='off');
            lgd.Layout.Tile = 'east';
            
            % Set plot width and plot height
            set(gcf, Units='centimeters');
            set(gcf, Position=[obj.margin obj.margin options.width options.height]);
        end

        function xy(obj,x,y,label,options)
            arguments (Input)
                obj (1,1) PlanarQuadrotorPlots
            end
            arguments (Input,Repeating)
                x (:,:) {mustBeVector,mustBeNumeric}
                y (:,:) {mustBeVector,mustBeNumeric}
                label (1,1) string {mustBeText}
            end
            arguments (Input)
                options.limits (:,:) cell {mustBeVector} = {''}
                options.colors (:,:) cell {mustBeVector} = {''}
                options.labels (:,:) cell {mustBeText,mustBeVector} = obj.xLabels
                options.labelFontSize (1,1) {mustBeNumeric} = obj.labelFontSize;
                options.tickFontSize (1,1) {mustBeNumeric,mustBeInteger} = obj.tickFontSize;
                options.legendFontSize (1,1) {mustBeNumeric,mustBeInteger} = obj.legendFontSize;
                options.lineWidth (1,1) {mustBeNumeric,mustBeNonnegative} = obj.lineWidth;
                options.alphaArea (1,1) {mustBeNumeric,mustBeNonnegative} = obj.alphaArea;
                options.width (1,1) {mustBeNumeric,mustBeNonnegative} = 12;
                options.height (1,1) {mustBeNumeric,mustBeNonnegative} = 10;
                options.limits_visible (1,1) string {mustBeMember(options.limits_visible,{'off','on'})} = "off"
            end
            
            % Get number of repeating arguments
            N = length(x);
            
            % Define colormaps (repeating)
            if any(cellfun(@isempty, options.colors))
                options.colors = obj.get_colormap(@cool,N,2);
            end
            
            % Create figure
            figure();
            tiledlayout('horizontal', TileSpacing='compact')
            nexttile; axis square, hold on; grid on;
            
            % Loop over repeating arguents
            for i = 1:N
                % Plot xy-trajectory
                plot(x{i},y{i}, ...
                     color=options.colors{i}{1}, ...
                     LineStyle='-', ...
                     LineWidth=options.lineWidth, ...
                     DisplayName=label{i});
            end
            
            % Plot limits
            if ~isempty(options.limits{1})
                xline(options.limits{1}, ...
                      color='k', ...
                      LineStyle='--', ...
                      LineWidth=1.5*options.lineWidth, ...
                      HandleVisibility='off', ...
                      Visible=options.limits_visible);
                xlim(options.limits{1})
            end
            if ~isempty(options.limits{2})
                yline(options.limits{2}, ...
                      color='k', ...
                      LineStyle='--', ...
                      LineWidth=1.5*options.lineWidth, ...
                      HandleVisibility='off', ...
                      Visible=options.limits_visible);
                ylim(options.limits{2})
            end
            
            % Plot legend
            lgd = legend(Interpreter='latex', ...
                         Orientation='vertical', ...
                         FontSize=options.legendFontSize, ...
                         Box='off');
            lgd.Layout.Tile = 'east';
            
            % Define labels and axis limits
            xlabel(['$',options.labels{1},'$'], Interpreter='latex', FontSize=options.labelFontSize);
            ylabel(['$',options.labels{2},'$'], Interpreter='latex', FontSize=options.labelFontSize);
            
            % Define fontsizes
            set(gca().XAxis, FontSize=options.tickFontSize)
            set(gca().YAxis, FontSize=options.tickFontSize)
            set(gca().XLabel, FontSize=options.labelFontSize)
            set(gca().YLabel, FontSize=options.labelFontSize)
            
            % Set plot width and plot height
            set(gcf, Units='centimeters');
            set(gcf, Position=[obj.margin obj.margin options.width options.height]);
        end

        function exportPlot(obj, fig, varargin)
            % Change figure and paper units
            set(fig, Units='centimeters');
            
            % Get width
            w_idx = find(cellfun(@(x) strcmp(x,'width'), varargin));
            if ~isempty(w_idx)
                width = varargin{w_idx+1};      %#ok
                varargin(w_idx:w_idx+1) = [];
            else
                pos = get(fig,'Position');
                width = pos(3);                 %#ok
            end
            
            % Get height
            h_idx = find(cellfun(@(x) strcmp(x,'height'), varargin));
            if ~isempty(h_idx)
                height = varargin{h_idx+1};
                varargin(h_idx:h_idx+1) = [];
            else
                pos = get(fig,'Position');
                height = pos(4);
            end
            
            % Set width and height
            set(fig, Position=[obj.margin,obj.margin,width,height]); %#ok
            
            % Change paper size (for export)
            set(fig, PaperUnits='centimeters');
            set(fig, PaperSize=[width height]);                      %#ok
            
            % Export plot
            saveas(fig, varargin{:})
        end

        function plot_states(obj,t,x,lower,upper,limits,colors,labels,lineWidth,lineStyle,alphaArea)
            % Plot only once
            for i=1:max(cellfun(@(x) size(x,1),x))
                % Plot dummy lines for legend
                plot(0, 0, 'k', ...
                     LineWidth=lineWidth, ...
                     LineStyle=lineStyle{i}, ...
                     DisplayName=['$',labels{i},'$'], ...
                     Visible='on');
                
                % Plot limits
                if ~isempty(limits{i})
                    yline(limits{i}, ...
                          color='k', ...
                          LineStyle=lineStyle{i}, ...
                          LineWidth=1.5*lineWidth, ...
                          HandleVisibility='off');
                end
            end
            
            % Loop over repeating arguments
            for i = 1:length(x)
                for j = 1:size(x{i},1)
                    % Plot tube
                    obj.tubes(t{i},lower{i}(j,:),upper{i}(j,:), ...
                              colors{i}{j}, ...
                              alphaArea, ...
                              ['$\mathcal{R}_{',labels{j},'}$'], ...
                              'off')
                    
                    % Plot trajectory
                    plot(t{i}, x{i}(j,:), ...
                         color=obj.darken(colors{i}{j},0.5), ...
                         LineWidth=lineWidth, ...
                         LineStyle=lineStyle{j}, ...
                         DisplayName=['$',labels{j},'$'], ...
                         HandleVisibility='off');
                end
            end
        end

        function plot_inputs(obj,t,u,lower,upper,limits,colors,labels,lineWidth,lineStyle,alphaArea)
            % Plot only once
            for i = 1:max(cellfun(@(x) size(x,1),u))
                % Plot dummy lines for legend
                plot(0, 0, 'k', ...
                     LineWidth=lineWidth, ...
                     LineStyle=lineStyle{i}, ...
                     DisplayName=['$',labels{i},'$'], ...
                     Visible='on');
                
                % Plot limits
                if ~isempty(limits{i})
                    yline(limits{i}, ...
                          color='k', ...
                          LineWidth=1.5*lineWidth, ...
                          LineStyle=lineStyle{i}, ...
                          HandleVisibility='off');
                end
            end
            
            % Loop over repeating arguments
            for i = 1:length(u)
                for j = 1:size(u{i},1)
                    % Plot tube
                    obj.discrete_tubes(t{i},lower{i}(j,:),upper{i}(j,:), ...
                                       colors{i}{j}, ...
                                       alphaArea, ...
                                       ['$\mathcal{R}_{',labels{j},'}$'], ...
                                       'off')
                    
                    % Plot trajectory
                    stairs(t{i}, [u{i}(j,1:end), u{i}(j,end)], ...
                           color=obj.darken(colors{i}{j},0.5), ...
                           LineWidth=lineWidth, ...
                           lineStyle=lineStyle{j}, ...
                           HandleVisibility='off');
                end
            end
        end
    end

    methods(Access=protected, Hidden=true, Static)
        function discrete_tubes(t,lower,upper,color,alphaArea,label,HandleVisibility)
            for j = 1:size(upper,2)
                x = [t(j),       t(j+1),   t(j+1),     t(j)];
                y = [lower(j), lower(j), upper(j), upper(j)];
                
                patch(x, y, color,...
                      FaceAlpha=alphaArea, ...
                      EdgeColor='none', ...
                      HandleVisibility=HandleVisibility, ...
                      DisplayName=label);
                HandleVisibility = 'off';
            end
        end

        function tubes(t,lower,upper,color,alphaArea,label,HandleVisibility)
            patch([t, fliplr(t)], [upper, fliplr(lower)], ...
                   color, ...
                   FaceAlpha=alphaArea, ...
                   EdgeColor='none', ...
                   HandleVisibility=HandleVisibility, ...
                   DisplayName=label);
        end

        function cmap = get_colormap(cFun,n,m)
            cmap = mat2cell(cFun(m*n),m*ones(1,n),3);
            if m > 1
                cmap = cellfun(@(X) num2cell(X,m), cmap, UniformOutput=false);
            end
        end

        function color = darken(color,beta)
            if beta > 0
                color = max([0,0,0], color - beta*[1,1,1]);
            else
                color = min([1,1,1], color - beta*[1,1,1]);
            end
        end
    end
end