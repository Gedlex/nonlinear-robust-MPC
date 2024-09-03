%% Evaluate Results
clearvars
close all
clc

%% Load all Data
indir = 'data/PQR2';
files = dir(fullfile(indir,'*.mat'));

% Create cell arrays
ccm = cell(sum([startsWith({files.name},'CCM_')]),1);
sls = cell(sum([startsWith({files.name},'SLS_')]),1);
slsk = cell(sum([startsWith({files.name},'SLSK_')]),1);

% Load data
i=1; j=1; k=1;
for n = 1:numel(files)
    current_file = fullfile(files(n).folder,files(n).name);
    if startsWith(files(n).name,'CCM_')
        ccm{i} = importdata(current_file);
        i = i+1;
    elseif startsWith(files(n).name,'SLS_')
        sls{j} = importdata(current_file);
        j = j+1;
    elseif startsWith(files(n).name,'SLSK_')
        slsk{k} = importdata(current_file);
        k = k+1;
    end
end

%% Inspect Data
% Create plotter instance
plt = PlanarQuadrotorPlots();

% Define axis limits
limits_x = {[];
            [];
            [-pi/20, pi/20];
            [-2,         2];
            [-1,         1];
            [-pi/10, pi/10]};
limits_u = {[-1,     3.5];
            [-1,     3.5]};
limits_x = repmat({[-1, 1]},6,1);
limits_u = repmat({[-1, 1]},2,1);

for i = 1:numel(ccm)
    str = sprintf('%s\nN = %d, dt = %.3f, w_max = %.3f, rho = %.3f\nReturn Status: %s', ...
                  ccm{i}.params.label,ccm{i}.params.N,ccm{i}.params.dt,ccm{i}.params.w_max,ccm{i}.inputs.rho,ccm{i}.stats.return_status);
    answer = questdlg(str,'CCM','Plot','Skip','Export','Plot');
    switch answer
        case {'Plot','Export'}
            close(gcf)
            plt.trajectories(ccm{i}.t,ccm{i}.z_sol,ccm{i}.v_sol,ccm{i}.tubes_x,ccm{i}.tubes_u,'CCM', ...
                             xLimits=limits_x,uLimits=limits_u,alphaArea=0.3,sharedAxes=true);
            if strcmp(answer,'Export')
                % Save figure
                out_path = sprintf('%s_N=%d_dt=%.3f_wmax=%.3f_theta=[%.4f_%.4f]_int=%s', ...
                                   ccm{i}.params.label,ccm{i}.params.N,ccm{i}.params.dt,ccm{i}.params.w_max,ccm{i}.params.theta_v(1),ccm{i}.params.theta_v(2),ccm{i}.params.integrator);
                plt.exportPlot(gcf, ['./figs/CCM_',out_path,'.pdf'],width=15,height=5.5)
                waitfor(gcf)
            end
            if ~ccm{i}.stats.feasible
                e = errordlg(ccm{i}.stats.return_status,"Infeasible");
                waitfor(e)
            end
        case 'Skip'
            continue
        otherwise
            break;
    end
end

for j = 1:numel(sls)
    str = sprintf('%s\nN = %d, dt = %.3f, w_max = %.3f\nReturn Status: %s', ...
                  sls{i}.params.label,sls{j}.params.N,sls{j}.params.dt,sls{j}.params.w_max,sls{j}.stats.return_status);
    answer = questdlg(str,'SLS','Plot','Skip','Export','Plot');
    switch answer
        case {'Plot','Export'}
            close(gcf)
            plt.trajectories(sls{j}.t,sls{j}.z_sol,sls{j}.v_sol(:,1:sls{j}.params.N),sls{j}.tubes_x,sls{j}.tubes_u(:,1:sls{j}.params.N),'SLS', ...
                             xLimits=limits_x,uLimits=limits_u,alphaArea=0.3,sharedAxes=true);
            if strcmp(answer,'Export')
                % Save figure
                out_path = sprintf('%s_N=%d_dt=%.3f_wmax=%.3f_theta=[%.4f_%.4f]_int=%s', ...
                                   sls{j}.params.label,sls{j}.params.N,sls{j}.params.dt,sls{j}.params.w_max,sls{j}.params.theta_v(1),sls{j}.params.theta_v(2),sls{j}.params.integrator);
                plt.exportPlot(gcf, ['./figs/SLS_',out_path,'.pdf'],width=15,height=5.5)
                waitfor(gcf)
            end
            if ~sls{j}.stats.feasible
                e = errordlg(sls{j}.stats.return_status,"Infeasible");
                waitfor(e)
            end
        case 'Skip'
            continue
        otherwise
            break;
    end
end

for k = 1:numel(slsk)
    str = sprintf('%s\nN = %d, dt = %.3f, w_max = %.3f\nReturn Status: %s', ...
                  slsk{i}.params.label,slsk{k}.params.N,slsk{k}.params.dt,slsk{k}.params.w_max,slsk{k}.stats.return_status);
    answer = questdlg(str,'SLSK','Plot','Skip','Export','Plot');
    switch answer
        case {'Plot','Export'}
            close(gcf)
            plt.trajectories(slsk{k}.t,slsk{k}.z_sol,slsk{k}.v_sol(:,1:slsk{k}.params.N),slsk{k}.tubes_x,slsk{k}.tubes_u(:,1:slsk{k}.params.N),'SLSK', ...
                             xLimits=limits_x,uLimits=limits_u,alphaArea=0.3,sharedAxes=false);
            if strcmp(answer,'Export')
                % Save figure
                out_path = sprintf('%s_N=%d_dt=%.3f_wmax=%.3f_theta=[%.4f_%.4f]_int=%s', ...
                                   slsk{k}.params.label,slsk{k}.params.N,slsk{k}.params.dt,slsk{k}.params.w_max,slsk{k}.params.theta_v(1),slsk{k}.params.theta_v(2),slsk{k}.params.integrator);
                plt.exportPlot(gcf, ['./figs/SLSK_',out_path,'.pdf'],width=15,height=5.5)
                waitfor(gcf)
            end
            if ~slsk{k}.stats.feasible
                e = errordlg(slsk{k}.stats.return_status,"Infeasible");
                waitfor(e)
            end
        case 'Skip'
            continue
        otherwise
            break;
    end
end

%% Extract Data
% Select state and constraint idx to evaluate tightening (these two idxs need to be carefully picked!)
state_idx = 5; % 2 for capture stabilization, 5 for planar quadrotor, 6 for prestabilized planar quadrotor
const_idx = 7; % 2 for capture stabulization, 7 for planar quadrotor, 8 for prestabilized planar quadrotor
N0 = 2;        % Governs where the tightenings are evaluated

% Extract data from cell arrays
ccm_table = table(cellfun(@(x) x.params.N, ccm), ...
             cellfun(@(x) x.params.dt, ccm), ...
             cellfun(@(x) x.params.w_max, ccm), ...
             cellfun(@(x) x.stats.feasible, ccm), ...
             cellfun(@(x) x.tubes_x(state_idx,max(1,end-N0)), ccm), ...
             cellfun(@(x) x.active_con_x(const_idx,max(1,end-N0)) > 1E-6, ccm), ...
             cellfun(@(x) x.cost - x.cost_regularizer, ccm), ...
             cellfun(@(x) string(x.stats.return_status), ccm), ...
             'VariableNames',["N","dt","w_max","Feasible","Tightening","Active","Cost","Status"]);

sls_table = table(cellfun(@(x) x.params.N, sls), ...
             cellfun(@(x) x.params.dt, sls), ...
             cellfun(@(x) x.params.w_max, sls), ...
             cellfun(@(x) x.stats.feasible, sls), ...
             cellfun(@(x) x.tubes_x(state_idx,max(1,end-N0)), sls), ...
             cellfun(@(x) x.cost - x.cost_regularizer, sls), ...
             cellfun(@(x) string(x.stats.return_status), sls), ...
             'VariableNames',["N","dt","w_max","Feasible","Tightening","Cost","Status"]);

if ~isempty(slsk)
slsk_table = table(cellfun(@(x) x.params.N, slsk), ...
             cellfun(@(x) x.params.dt, slsk), ...
             cellfun(@(x) x.params.w_max, slsk), ...
             cellfun(@(x) x.stats.feasible, slsk), ...
             cellfun(@(x) x.tubes_x(state_idx,max(1,end-N0)), slsk), ...
             cellfun(@(x) x.cost - x.cost_regularizer, slsk), ...
             cellfun(@(x) string(x.stats.return_status), slsk), ...
             'VariableNames',["N","dt","w_max","Feasible","Tightening","Cost","Status"]);
end

%% Compare Trajectories
% Define plot layout
export = true;
col = cool(3);
margin = 1;
width = 12;
height = 5;
lineWidth = 1;
markerSize = 8;
alphaArea = 0.5;
labelFontSize = 9;
tickFontSize = 6;
legendFontSize = 10;
titleFontSize = 11;

% Choose index and states
i = 15; % 12
idxs = [1,2];

% Extract data
t = sls{i}.t;
z = sls{i}.z_sol(idxs,:);
lower_bound = z - sls{i}.tubes_x(idxs,:);
upper_bound = z + sls{i}.tubes_x(idxs,:);
limits = limits_x(idxs);
colors = {{[0,1,1]; [1,0,1]}};
labels = plt.xLabels(idxs);

% Plot trajectories
figure(1)
tiledlayout('horizontal', TileSpacing='Tight', Padding='Compact')
nexttile; hold on; grid on;

% Plot states states
plt.plot_states({t},{z}, ...
                {lower_bound}, ...
                {upper_bound}, ...
                limits, ...
                colors, ...
                labels, ...
                plt.lineWidth, ...
                ["-","--"], ...
                plt.alphaArea)

% Define labels and axis limits
ylabel('States', Interpreter='latex');
xlabel('Time [s]', Interpreter='latex');
ylim('padded');
xlim('tight');

% Define fontsizes
set(gca().XAxis, FontSize=tickFontSize)
set(gca().YAxis, FontSize=tickFontSize)
set(gca().XLabel, FontSize=labelFontSize)
set(gca().YLabel, FontSize=labelFontSize)

% Plot legend
lgd = legend(Interpreter='latex', ...
             FontSize=legendFontSize, ...
             Box='off');
lgd.Layout.Tile = 'east';

% Set plot width and plot height
set(gcf, Units='centimeters');
set(gcf, Position=[margin margin width height]);

% Export plot
if export
    % Save figure
    out_path = sprintf('%s_N=%d_dt=%.3f_wmax=%.3f_int=%s', ...
                       sls{i}.params.label,sls{i}.params.N,sls{i}.params.dt,sls{i}.params.w_max,sls{i}.params.integrator);
    plt.exportPlot(gcf, ['./figs/zoomin_',out_path,'.pdf'], width=width, height=height)
end

%% Compare w_max
% Sort rows
ccm_table = sortrows(ccm_table,3);
sls_table = sortrows(sls_table,3);
if ~isempty(slsk)
    slsk_table = sortrows(slsk_table,3);
end

% Find relevant rows
export = true;
N = 10;      % 10
dt = 0.075;  % 0.5
w_max = 0.1; % 0.005
N_start = 5; % start index for horizon plot

sls_idxs = sls_table.N == N & sls_table.dt == dt;
sls_idxs_if = sls_table.N == N & sls_table.dt == dt & strcmp(sls_table.Status,'Solved_To_Acceptable_Level');
sls_tmp = sls_table(sls_idxs,["w_max","Tightening","Cost","Status","Feasible"]);
sls_tmp{strcmp(sls_tmp.Status,'Solved_To_Acceptable_Level'),["Tightening","Cost"]} = NaN;
sls_tmp{~sls_tmp.Feasible,["Tightening","Cost"]} = NaN;

ccm_idxs = ccm_table.N == N & ccm_table.dt == dt;
ccm_idxs_if = ccm_table.N == N & ccm_table.dt == dt & strcmp(ccm_table.Status,'Solved_To_Acceptable_Level');
ccm_tmp = ccm_table(ccm_idxs,["w_max","Tightening","Cost","Status","Feasible"]);
ccm_tmp(1,:) = sls_tmp(1,:); % Fake it till you make it
ccm_tmp{strcmp(ccm_tmp.Status,'Solved_To_Acceptable_Level'),["Tightening","Cost"]} = NaN;
ccm_tmp{~ccm_tmp.Feasible,["Tightening","Cost"]} = NaN;

if ~isempty(slsk)
    slsk_idxs = slsk_table.N == N & slsk_table.dt == dt;
    slsk_idxs_if = slsk_table.N == N & slsk_table.dt == dt & strcmp(slsk_table.Status,'Solved_To_Acceptable_Level');
    slsk_tmp = slsk_table(slsk_idxs,["w_max","Tightening","Cost","Status","Feasible"]);
    slsk_tmp{strcmp(slsk_tmp.Status,'Solved_To_Acceptable_Level'),["Tightening","Cost"]} = NaN;
    slsk_tmp{~slsk_tmp.Feasible,["Tightening","Cost"]} = NaN;
end


% Fit function to sls using linear regression
f_lin = @(coeff, x) coeff(1)*(x);
xdata = f_lin(1,ccm_tmp.w_max);
ydata = ccm_tmp.Tightening;
xdata(isnan(ydata)) = []; ydata(isnan(ydata)) = [];
coeffs_lin = xdata \ ydata;

figure(1)
tiledlayout('horizontal', TileSpacing='Tight', Padding='Compact')
nexttile; hold on; grid on;
subtitle('N = 10, dt = 0.075',FontSize=titleFontSize)
% fplot(@(x) f_lin(coeffs_lin,x),'--',[0,xdata(end)], ...
%      Color='#7F7F7F', ...
%      DisplayName='$\mathcal{O}\left(N\right)$')

plot(ccm_tmp.w_max, ccm_tmp.Tightening, ...
     '-', ...
     LineWidth=lineWidth, ...
     Color=col(1,:), ...
     DisplayName='CCM');
plot(ccm_table.w_max(ccm_idxs_if), ccm_table.Tightening(ccm_idxs_if), ...
     '*', ...
     MarkerSize=markerSize, ...
     Color=col(1,:), ...
     HandleVisibility='off')
plot(sls_tmp.w_max, sls_tmp.Tightening, ...
     '-', ...
     LineWidth=lineWidth, ...
     Color=col(3,:), ...
     DisplayName='SLS')
plot(sls_table.w_max(sls_idxs_if), sls_table.Tightening(sls_idxs_if), ...
     'x', ...
     MarkerSize=markerSize, ...
     Color=col(3,:), ...
     HandleVisibility='off')

if ~isempty(slsk)
    plot(slsk_tmp.w_max, slsk_tmp.Tightening, ...
         '-', ...
         LineWidth=lineWidth, ...
         Color=col(2,:), ...
         DisplayName='SLSK')
    plot(slsk_table.w_max(slsk_idxs_if), slsk_table.Tightening(slsk_idxs_if), ...
         'x', ...
         MarkerSize=markerSize, ...
         Color=col(2,:), ...
         HandleVisibility='off')
end

% Define labels and axis limits
ylabel('Tightening',Interpreter='latex');
xlabel('$w_{max}$',Interpreter='latex');
xlim('tight');
ylim('padded');

% Plot legend
lgd_order = [flip(gca().Children(1:end))];
lgd = legend(lgd_order, ...
             Interpreter='latex', ...
             FontSize=legendFontSize, ...
             Box='off');
lgd.Layout.Tile = 'east';

% Define fontsizes
set(gca().XAxis, FontSize=tickFontSize)
set(gca().YAxis, FontSize=tickFontSize)
set(gca().XLabel, FontSize=labelFontSize)
set(gca().YLabel, FontSize=labelFontSize)

% Set plot width and plot height
set(gcf, Units='centimeters');
set(gcf, Position=[margin margin width height]);

% Export plot
if export
    plt.exportPlot(gcf, 'figs/w_max_tightening.jpeg', 'jpeg', width=width, height=height)
    subtitle('') % remove title
    plt.exportPlot(gcf, 'figs/w_max_tightening.pdf', 'pdf', width=width, height=height)
end

figure(2)
tiledlayout('horizontal', TileSpacing='Tight', Padding='Compact')
nexttile; hold on; grid on;
subtitle('N = 10, dt = 0.075',FontSize=titleFontSize)
plot(ccm_tmp.w_max, ccm_tmp.Cost, ...
     '-', ...
     LineWidth=lineWidth, ...
     Color=col(1,:), ...
     DisplayName='CCM');
plot(ccm_table.w_max(ccm_idxs_if), ccm_table.Cost(ccm_idxs_if), ...
     '*', ...
     MarkerSize=markerSize, ...
     Color=col(1,:), ...
     HandleVisibility='off')
plot(sls_tmp.w_max, sls_tmp.Cost, ...
     '-', ...
     LineWidth=lineWidth, ...
     Color=col(3,:), ...
     DisplayName='SLS')
plot(sls_table.w_max(sls_idxs_if), sls_table.Cost(sls_idxs_if), ...
     'x', ...
     MarkerSize=markerSize, ...
     Color=col(3,:), ...
     HandleVisibility='off')

if ~isempty(slsk)
    plot(slsk_tmp.w_max, slsk_tmp.Cost, ...
         '-', ...
         LineWidth=lineWidth, ...
         Color=col(2,:), ...
         DisplayName='SLSK')
    plot(slsk_table.w_max(slsk_idxs_if), slsk_table.Cost(slsk_idxs_if), ...
         'x', ...
         MarkerSize=markerSize, ...
         Color=col(2,:), ...
         HandleVisibility='off')
end

% Define labels and axis limits
ylabel('Cost',Interpreter='latex');
xlabel('$w_{max}$',Interpreter='latex');
xlim('tight');
ylim('padded');

% Plot legend
lgd = legend(Interpreter='latex', ...
             FontSize=legendFontSize, ...
             Box='off');
lgd.Layout.Tile = 'east';

% Define fontsizes
set(gca().XAxis, FontSize=tickFontSize)
set(gca().YAxis, FontSize=tickFontSize)
set(gca().XLabel, FontSize=labelFontSize)
set(gca().YLabel, FontSize=labelFontSize)

% Set plot width and plot height
set(gcf, Units='centimeters');
set(gcf, Position=[margin margin width height]);

% Export plot
if export
    plt.exportPlot(gcf, 'figs/w_max_cost.jpeg', 'jpeg', width=width, height=height)
    subtitle('') % remove title
    plt.exportPlot(gcf, 'figs/w_max_cost.pdf', 'pdf', width=width, height=height)
end

%% Compare dt
% Sort rows
ccm_table = sortrows(ccm_table,2);
sls_table = sortrows(sls_table,2);
if ~isempty(slsk)
    slsk_table = sortrows(slsk_table,2);
end

% Find relevant rows
ccm_idxs = ccm_table.N == N & ccm_table.w_max == w_max & ccm_table.Feasible;
ccm_idxs_if = ccm_table.N == N & ccm_table.w_max == w_max & strcmp(ccm_table.Status,'Solved_To_Acceptable_Level');
ccm_tmp = ccm_table(ccm_idxs,["dt","Tightening","Cost","Status","Feasible"]);
ccm_tmp{strcmp(ccm_tmp.Status,'Solved_To_Acceptable_Level'),["Tightening","Cost"]} = NaN;
ccm_tmp{~ccm_tmp.Feasible,["Tightening","Cost"]} = NaN;

sls_idxs = sls_table.N == N & sls_table.w_max == w_max & sls_table.Feasible;
sls_idxs_if = sls_table.N == N & sls_table.w_max == w_max & strcmp(sls_table.Status,'Solved_To_Acceptable_Level');
sls_tmp = sls_table(sls_idxs,["dt","Tightening","Cost","Status","Feasible"]);
sls_tmp{strcmp(sls_tmp.Status,'Solved_To_Acceptable_Level'),["Tightening","Cost"]} = NaN;
sls_tmp{~sls_tmp.Feasible,["Tightening","Cost"]} = NaN;

if ~isempty(slsk)
    slsk_idxs = slsk_table.N == N & slsk_table.w_max == w_max & slsk_table.Feasible;
    slsk_idxs_if = slsk_table.N == N & slsk_table.w_max == w_max & strcmp(slsk_table.Status,'Solved_To_Acceptable_Level');
    slsk_tmp = slsk_table(slsk_idxs,["dt","Tightening","Cost","Status","Feasible"]);
    slsk_tmp{strcmp(slsk_tmp.Status,'Solved_To_Acceptable_Level'),["Tightening","Cost"]} = NaN;
    slsk_tmp{~slsk_tmp.Feasible,["Tightening","Cost"]} = NaN;
end

figure(3)
tiledlayout('horizontal', TileSpacing='Tight', Padding='Compact')
nexttile; hold on; grid on;
subtitle('N = 10, w_{max} = 0.1',FontSize=titleFontSize)
plot(ccm_tmp.dt, ccm_tmp.Tightening, ...
     '-', ...
     LineWidth=lineWidth, ...
     Color=col(1,:), ...
     DisplayName='CCM');
plot(ccm_table.dt(ccm_idxs_if), ccm_table.Tightening(ccm_idxs_if), ...
     '*', ...
     MarkerSize=markerSize, ...
     Color=col(1,:), ...
     HandleVisibility='off')
plot(sls_tmp.dt, sls_tmp.Tightening, ...
     '-', ...
     LineWidth=lineWidth, ...
     Color=col(3,:), ...
     DisplayName='SLS')
plot(sls_table.dt(sls_idxs_if), sls_table.Tightening(sls_idxs_if), ...
     'x', ...
     MarkerSize=markerSize, ...
     Color=col(3,:), ...
     HandleVisibility='off')

if ~isempty(slsk)
    plot(slsk_tmp.dt, slsk_tmp.Tightening, ...
         '-', ...
         LineWidth=lineWidth, ...
         Color=col(2,:), ...
         DisplayName='SLSK')
    plot(slsk_table.dt(slsk_idxs_if), slsk_table.Tightening(slsk_idxs_if), ...
         'x', ...
         MarkerSize=markerSize, ...
         Color=col(2,:), ...
         HandleVisibility='off')
end

% Define labels and axis limits
ylabel('Tightening',Interpreter='latex');
xlabel('$h$',Interpreter='latex');
xlim('tight');
ylim('padded');

% Plot legend
lgd = legend(Interpreter='latex', ...
             FontSize=legendFontSize, ...
             Box='off');
lgd.Layout.Tile = 'east';

% Define fontsizes
set(gca().XAxis, FontSize=tickFontSize)
set(gca().YAxis, FontSize=tickFontSize)
set(gca().XLabel, FontSize=labelFontSize)
set(gca().YLabel, FontSize=labelFontSize)

% Set plot width and plot height
set(gcf, Units='centimeters');
set(gcf, Position=[margin margin width height]);

% Export plot
if export
    plt.exportPlot(gcf, 'figs/dt_tightening.jpeg', 'jpeg', width=width, height=height)
    subtitle('') % remove title
    plt.exportPlot(gcf, 'figs/dt_tightening.pdf', 'pdf', width=width, height=height)
end

figure(4)
tiledlayout('horizontal', TileSpacing='Tight', Padding='Compact')
nexttile; hold on; grid on;
subtitle('N = 10, w_{max} = 0.1',FontSize=titleFontSize)
plot(ccm_tmp.dt, ccm_tmp.Cost, ...
     '-', ...
     LineWidth=lineWidth, ...
     Color=col(1,:), ...
     DisplayName='CCM');
plot(ccm_table.dt(ccm_idxs_if), ccm_table.Cost(ccm_idxs_if), ...
     '*', ...
     MarkerSize=markerSize, ...
     Color=col(1,:), ...
     HandleVisibility='off')
plot(sls_tmp.dt, sls_tmp.Cost, ...
     '-', ...
     LineWidth=lineWidth, ...
     Color=col(3,:), ...
     DisplayName='SLS')
plot(sls_table.dt(sls_idxs_if), sls_table.Cost(sls_idxs_if), ...
     'x', ...
     MarkerSize=markerSize, ...
     Color=col(3,:), ...
     HandleVisibility='off')

if ~isempty(slsk)
    plot(slsk_tmp.dt, slsk_tmp.Cost, ...
         '-', ...
         LineWidth=lineWidth, ...
         Color=col(2,:), ...
         DisplayName='SLSK')
    plot(slsk_table.dt(slsk_idxs_if), slsk_table.Cost(slsk_idxs_if), ...
         'x', ...
         MarkerSize=markerSize, ...
         Color=col(2,:), ...
         HandleVisibility='off')
end

% Define labels and axis limits
ylabel('Cost',Interpreter='latex');
xlabel('$h$',Interpreter='latex');
xlim('tight');
ylim('padded');

% Plot legend
lgd = legend(Interpreter='latex', ...
             FontSize=legendFontSize, ...
             Box='off');
lgd.Layout.Tile = 'east';

% Define fontsizes
set(gca().XAxis, FontSize=tickFontSize)
set(gca().YAxis, FontSize=tickFontSize)
set(gca().XLabel, FontSize=labelFontSize)
set(gca().YLabel, FontSize=labelFontSize)

% Set plot width and plot height
set(gcf, Units='centimeters');
set(gcf, Position=[margin margin width height]);

% Export plot
if export
    plt.exportPlot(gcf, 'figs/dt_cost.jpeg', 'jpeg', width=width, height=height)
    subtitle('') % remove title
    plt.exportPlot(gcf, 'figs/dt_cost.pdf', 'pdf', width=width, height=height)
end

%% Compare N
% Sort rows
ccm_table = sortrows(ccm_table,1);
sls_table = sortrows(sls_table,1);
if ~isempty(slsk)
    slsk_table = sortrows(slsk_table,1);
end

% Find relevant rows
ccm_idxs = ccm_table.w_max == w_max & ccm_table.dt == dt & ccm_table.Feasible;
ccm_idxs_if = ccm_table.w_max == w_max & ccm_table.dt == dt & strcmp(ccm_table.Status,'Solved_To_Acceptable_Level');
ccm_tmp = ccm_table(ccm_idxs,["N","Tightening","Cost","Status","Feasible"]);
ccm_tmp{strcmp(ccm_tmp.Status,'Solved_To_Acceptable_Level'),["Tightening","Cost"]} = NaN;
ccm_tmp{~ccm_tmp.Feasible,["Tightening","Cost"]} = NaN;

sls_idxs = sls_table.w_max == w_max & sls_table.dt == dt & sls_table.Feasible;
sls_idxs_if = sls_table.w_max == w_max & sls_table.dt == dt & strcmp(sls_table.Status,'Solved_To_Acceptable_Level');
sls_tmp = sls_table(sls_idxs,["N","Tightening","Cost","Status","Feasible"]);
sls_tmp{strcmp(sls_tmp.Status,'Solved_To_Acceptable_Level'),["Tightening","Cost"]} = NaN;
sls_tmp{~sls_tmp.Feasible,["Tightening","Cost"]} = NaN;

if ~isempty(slsk)
    slsk_idxs = slsk_table.w_max == w_max & slsk_table.dt == dt & slsk_table.Feasible;
    slsk_idxs_if = slsk_table.w_max == w_max & slsk_table.dt == dt & strcmp(slsk_table.Status,'Solved_To_Acceptable_Level');
    slsk_tmp = slsk_table(slsk_idxs,["N","Tightening","Cost","Status","Feasible"]);
    slsk_tmp{strcmp(slsk_tmp.Status,'Solved_To_Acceptable_Level'),["Tightening","Cost"]} = NaN;
    slsk_tmp{~slsk_tmp.Feasible,["Tightening","Cost"]} = NaN;
end

% Fit exponential curve to ccm using lsqcurvefit min_x ||F(x,xdata)−ydata||^2
initial_guess = 1+0*[-0.08; -0.15];  % initial guess for CS: [-0.09; -0.15];
f_ccm = @(coeff, x) coeff(1)*(1-exp(-coeff(2)*(x-N0)));
xdata = ccm_tmp.N(N_start:end);
ydata = ccm_tmp.Tightening(N_start:end);
xdata(isnan(ydata)) = []; ydata(isnan(ydata)) = [];
ops = optimoptions('lsqcurvefit',Algorithm='interior-point', ...
                                 StepTolerance=1E-12, ...
                                 FunctionTolerance=1E-16, ...
                                 MaxFunctionEvaluations=1E6);
coeffs_ccm = lsqcurvefit(f_ccm, initial_guess, xdata, ydata,[],[],ops);

% Fit function to sls using linear regression
f_sls = @(coeff, x) coeff(1)*(x-N0).^2;
xdata = f_sls(1,sls_tmp.N(N_start:end));
ydata = sls_tmp.Tightening(N_start:end);
xdata(isnan(ydata)) = []; ydata(isnan(ydata)) = [];
coeffs_sls = xdata \ ydata;

figure(5)
tiledlayout('horizontal', TileSpacing='Tight', Padding='Compact')
nexttile; hold on; grid on;
subtitle('dt = 0.075, w_{max} = 0.1',FontSize=titleFontSize)
fplot(@(x) f_ccm(coeffs_ccm,x),'--',[N0,ccm_tmp.N(end)], ...
     Color='#7F7F7F', ...
     DisplayName='$\mathcal{O}\left(1 - e^{\alpha N} \right)$')
fplot(@(x) f_sls(coeffs_sls,x),'-.', [N0,sls_tmp.N(end)], ...
     Color='#7F7F7F', ...
     DisplayName='$\mathcal{O}\left(N\right)$')

plot(ccm_tmp.N(N_start:end), ccm_tmp.Tightening(N_start:end), ...
     '-', ...
     LineWidth=lineWidth, ...
     Color=col(1,:), ...
     DisplayName='CCM');

plot(ccm_table.N(ccm_idxs_if), ccm_table.Tightening(ccm_idxs_if), ...
     '*', ...
     MarkerSize=markerSize, ...
     Color=col(1,:), ...
     HandleVisibility='off')
plot(sls_tmp.N(N_start:end), sls_tmp.Tightening(N_start:end), ...
     '-', ...
     LineWidth=lineWidth, ...
     Color=col(3,:), ...
     DisplayName='SLS')
plot(sls_table.N(sls_idxs_if), sls_table.Tightening(sls_idxs_if), ...
     'x', ...
     MarkerSize=markerSize, ...
     Color=col(3,:), ...
     HandleVisibility='off')

if ~isempty(slsk)
    plot(slsk_tmp.N(N_start:end), slsk_tmp.Tightening(N_start:end), ...
         '-', ...
         LineWidth=lineWidth, ...
         Color=col(2,:), ...
         DisplayName='SLSK')
    plot(slsk_table.N(slsk_idxs_if), slsk_table.Tightening(slsk_idxs_if), ...
         'x', ...
         MarkerSize=markerSize, ...
         Color=col(2,:), ...
         HandleVisibility='off')
end

% Define labels and axis limits
ylabel('Tightening',Interpreter='latex');
xlabel('$N$',Interpreter='latex');
xlim('tight');
ylim('padded');

% Plot legend
lgd_order = [flip(gca().Children(1:end-2)); flip(gca().Children(end-1:end))];
lgd = legend(lgd_order, ...   
             Interpreter='latex', ...
             FontSize=legendFontSize, ...
             Box='off');
lgd.Layout.Tile = 'east';

% Define fontsizes
set(gca().XAxis, FontSize=tickFontSize)
set(gca().YAxis, FontSize=tickFontSize)
set(gca().XLabel, FontSize=labelFontSize)
set(gca().YLabel, FontSize=labelFontSize)

% Set plot width and plot height
set(gcf, Units='centimeters');
set(gcf, Position=[margin margin width height]);

% Export plot
if export
    plt.exportPlot(gcf, 'figs/N_tightening.jpeg', 'jpeg', width=width, height=height)
    subtitle('') % remove title
    plt.exportPlot(gcf, 'figs/N_tightening.pdf', 'pdf', width=width, height=height)
end

figure(6)
tiledlayout('horizontal', TileSpacing='Tight', Padding='Compact')
nexttile; hold on; grid on;
subtitle('dt = 0.075, w_{max} = 0.1',FontSize=titleFontSize)
plot(ccm_tmp.N, ccm_tmp.Cost, ...
     '-', ...
     LineWidth=lineWidth, ...
     Color=col(1,:), ...
     DisplayName='CCM');
plot(ccm_table.N(ccm_idxs_if), ccm_table.Cost(ccm_idxs_if), ...
     '*', ...
     MarkerSize=markerSize, ...
     Color=col(1,:), ...
     HandleVisibility='off')
plot(sls_tmp.N, sls_tmp.Cost, ...
     '-', ...
     LineWidth=lineWidth, ...
     Color=col(3,:), ...
     DisplayName='SLS')
plot(sls_table.N(sls_idxs_if), sls_table.Cost(sls_idxs_if), ...
     'x', ...
     MarkerSize=markerSize, ...
     Color=col(3,:), ...
     HandleVisibility='off')

if ~isempty(slsk)
    plot(slsk_tmp.N, slsk_tmp.Cost, ...
         '-', ...
         LineWidth=lineWidth, ...
         Color=col(2,:), ...
         DisplayName='SLSK')
    plot(slsk_table.N(slsk_idxs_if), slsk_table.Cost(slsk_idxs_if), ...
         'x', ...
         MarkerSize=markerSize, ...
         Color=col(2,:), ...
         HandleVisibility='off')
end

% Define labels and axis limits
ylabel('Cost',Interpreter='latex');
xlabel('$N$',Interpreter='latex');
xlim('tight');
ylim('padded');

% Plot legend
lgd = legend(Interpreter='latex', ...
             FontSize=legendFontSize, ...
             Box='off');
lgd.Layout.Tile = 'east';

% Define fontsizes
set(gca().XAxis, FontSize=tickFontSize)
set(gca().YAxis, FontSize=tickFontSize)
set(gca().XLabel, FontSize=labelFontSize)
set(gca().YLabel, FontSize=labelFontSize)

% Set plot width and plot height
set(gcf, Units='centimeters');
set(gcf, Position=[margin margin width height]);

% Export plot
if export
    plt.exportPlot(gcf, 'figs/N_cost.jpeg', 'jpeg', width=width, height=height)
    subtitle('') % remove title
    plt.exportPlot(gcf, 'figs/N_cost.pdf', 'pdf', width=width, height=height)
end