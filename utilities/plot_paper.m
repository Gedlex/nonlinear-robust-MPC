% -----------------------------------------------------------------------------
% File: plot_paper.m
% Author: Antoine Leeman (aleeman@ethz.ch)
% Date: 09th September 2023
% License: MIT
% Reference:
%{
@inproceedings{leeman2023_CDC,
title={Robust Optimal Control for Nonlinear Systems with Parametric Uncertainties via System Level Synthesis},
author={Leeman, Antoine P. and Sieber, Jerome and Bennani, Samir and Zeilinger, Melanie N.},
booktitle = {Proc. of the 62nd IEEE Conf. on Decision and Control (CDC)},
doi={10.1109/CDC49753.2023.10383271},
pages={4784-4791},
year={2023}}
%}
% Link: https://arxiv.org/abs/2304.00752
% -----------------------------------------------------------------------------

%% Compute Trajectories
% Compute Nominal MPC (used for comparison)
[z_init,v_init] = sls.computeNominalMPC();

% Get feedback matrix
K = M_sol/R_sol;
theta_true = params.theta_true;

% Predict closed loop trajectory with open loop control strategy
x_cl = zeros(sls.nx,sls.N+1);
u_cl = zeros(sls.nu,sls.N+1);
x_cl(:,1) = params.x0;
for k = 1 : sls.N
    Delta_x = reshape(x_cl(:,2:end) -s.z_sol(:,2:end), [sls.N*sls.nx,1]);
    u_cl = reshape([zeros(sls.nu,1);K*Delta_x],[sls.nu,sls.N+1]) + s.v_sol;
    x_cl(:,k+1)= sys.ddyn(x_cl(:,k),u_cl(:,k),theta=theta_true,d=params.w_max*(2*(rand(params.nw,1)>0.5)-1));
end

% Predict nominal trajectory with open loop control strategy
x_cl_init = zeros(sls.nx,sls.N+1);
u_cl_init = zeros(sls.nu,sls.N+1);
x_cl_init(:,1) = params.x0;
for k = 1:sls.N
    u_cl_init(:,k) = v_init(:,k);
    x_cl_init(:,k+1)= sys.ddyn(x_cl_init(:,k),u_cl_init(:,k),theta=theta_true,d=params.w_max*(2*(rand(params.nw,1)>0.5)-1));
end

%% Figure 1
figure(1);
clf
colormap = [0.0504    0.0298    0.5280
            0.4934    0.0115    0.6580
            0.7964    0.2780    0.4713
            0.9722    0.5817    0.2541
            0.9400    0.9752    0.1313];

color1 = colormap(1,:);
color2 = colormap(3,:);
color3 = colormap(2,:);
color4 = colormap(4,:);

labelFontSize = 6;
tickFontSize = 8;
legendfontsize = 10;
alpha_line = 0.9;
alpha_area = 0.5;

highlightColor = [hex2dec('7F') hex2dec('00') hex2dec('FF')] / 255;
highlightAlpha = 0.2;

subplot(1,4,1);
hold on;

upper_y = z_sol(1,:)+tubes_x(1,:);
lower_y = z_sol(1,:)-tubes_x(1,:);
h_rx = fill([0:params.dt:params.T, fliplr(0:params.dt:params.T)], [upper_y, fliplr(lower_y)],[0.8, 0.8, 0.8], 'FaceColor',color1, 'EdgeColor', 'none', 'FaceAlpha', alpha_area);

upper_y = z_sol(2,:)+tubes_x(2,:);
lower_y = z_sol(2,:)-tubes_x(2,:);
h_ry = fill([0:params.dt:params.T, fliplr(0:params.dt:params.T)], [upper_y, fliplr(lower_y)], [0.8, 0.8, 0.8],'FaceColor',color2, 'EdgeColor', 'none', 'FaceAlpha', alpha_area);

plot(0:params.dt:params.T, x_cl(2,:),'k--', 'LineWidth', 1.1, 'MarkerSize', 10);
plot(0:params.dt:params.T, x_cl(1,:),'k--', 'LineWidth', 1.1, 'MarkerSize', 10);
plot(0:params.dt:params.T, x_cl_init(2,:),'k:', 'LineWidth', 1.1, 'MarkerSize', 10);
plot(0:params.dt:params.T, x_cl_init(1,:),'k:', 'LineWidth', 1.1, 'MarkerSize', 10);

h_x = plot(0:params.dt:params.T, z_sol(1,:),'color',color1,'LineStyle','-', 'LineWidth', 1.5, 'MarkerSize', 10);
h_y = plot(0:params.dt:params.T, z_sol(2,:),'color',color2,'LineStyle','-', 'LineWidth', 1.5, 'MarkerSize', 10);

yline(1, '-k', 'LineWidth', 1.5);
yline(-1, '-k', 'LineWidth', 1.5);
lgd1 = legend([h_x, h_y, h_rx, h_ry], {'$p_x^\star$', '$p_y^\star$','$\mathcal{R}_{p_x^\star}$','$\mathcal{R}_{p_y^\star}$'}, 'Interpreter', 'latex', 'Location', 'best','fontsize',legendfontsize,'Position',[0.1712 0.1888 0.0763 0.2338]);
set(lgd1, 'Box', 'off', 'Color', 'none');

xlim([0, 5]);
ylim([-1, 1.1]);
ylabel('states', 'FontSize', labelFontSize);
xlabel('time [-]', 'FontSize', labelFontSize);
set(gca, 'FontSize', tickFontSize);
grid on;

subplot(1,4,2);
hold on;

upper_y = z_sol(3,:)+tubes_x(3,:);
lower_y = z_sol(3,:)-tubes_x(3,:);
rx = fill([0:params.dt:params.T, fliplr(0:params.dt:params.T)], [upper_y, fliplr(lower_y)],[0.8, 0.8, 0.8], 'FaceColor',color1, 'EdgeColor', 'none', 'FaceAlpha', 0.5);

upper_y = z_sol(6,:)+tubes_x(6,:);
lower_y = z_sol(6,:)-tubes_x(6,:);
ry = fill([0:params.dt:params.T, fliplr(0:params.dt:params.T)], [upper_y, fliplr(lower_y)],[0.8, 0.8, 0.8], 'FaceColor',color2, 'EdgeColor', 'none', 'FaceAlpha', 0.5);


h_x = plot(0:params.dt:params.T, z_sol(3,:),'color',color1,'LineStyle','-', 'LineWidth', 1.5, 'MarkerSize', 10);
h_y = plot(0:params.dt:params.T, z_sol(6,:),'color',color2,'LineStyle','-', 'LineWidth', 1.5, 'MarkerSize', 10);
plot(0:params.dt:params.T, x_cl(3,:),'k--','LineWidth', 1.1, 'MarkerSize', 10);
plot(0:params.dt:params.T, x_cl(6,:),'k--', 'LineWidth', 1.1, 'MarkerSize', 10);
plot(0:params.dt:params.T, x_cl_init(3,:),'k:','LineWidth', 1.1, 'MarkerSize', 10);
plot(0:params.dt:params.T, x_cl_init(6,:),'k:', 'LineWidth', 1.1, 'MarkerSize', 10);

yline(1, '-k', 'LineWidth', 1.5);
yline(-1, '-k', 'LineWidth', 1.5);
lgd2 = legend([h_x, h_y, rx, ry], {'$\psi^\star$', '$\dot \psi^\star$','$\mathcal{R}_{\psi^\star}$','$\mathcal{R}_{\dot \psi^\star}$'}, 'Interpreter', 'latex', 'Location', 'best','fontsize',legendfontsize,'Position',[0.3925 0.6182 0.0777 0.2379]);
set(lgd2, 'Box', 'off', 'Color', 'none');
xlim([0, 5]);
ylim([-1, 1.1]);
xlabel('time [-]', 'FontSize', labelFontSize);
set(gca,'YTickLabel', [], 'FontSize', tickFontSize);

grid on;
subplot(1,4,3);

hold on
grid on;
upper_y = z_sol(4,:)+tubes_x(4,:);
lower_y = z_sol(4,:)-tubes_x(4,:);
rx = fill([0:params.dt:params.T, fliplr(0:params.dt:params.T)], [upper_y, fliplr(lower_y)],[0.8, 0.8, 0.8], 'FaceColor',color1, 'EdgeColor', 'none', 'FaceAlpha', 0.5);

upper_y = z_sol(5,:)+tubes_x(5,:);
lower_y = z_sol(5,:)-tubes_x(5,:);
ry = fill([0:params.dt:params.T, fliplr(0:params.dt:params.T)], [upper_y, fliplr(lower_y)],[0.8, 0.8, 0.8], 'FaceColor',color2, 'EdgeColor', 'none', 'FaceAlpha', 0.5);

plot(0:params.dt:params.T, x_cl(4,:),'k--', 'LineWidth', 1.5, 'MarkerSize', 10);
plot(0:params.dt:params.T, x_cl(5,:),'k--','LineWidth', 1.5, 'MarkerSize', 10);
plot(0:params.dt:params.T, x_cl_init(4,:),'k:', 'LineWidth', 1.5, 'MarkerSize', 10);
plot(0:params.dt:params.T, x_cl_init(5,:),'k:','LineWidth', 1.5, 'MarkerSize', 10);

h_x = plot(0:params.dt:params.T, z_sol(4,:),'color',color1,'LineStyle','-', 'LineWidth', 1.1, 'MarkerSize', 10);
h_y = plot(0:params.dt:params.T, z_sol(5,:),'color',color2,'LineStyle','-', 'LineWidth', 1.1, 'MarkerSize', 10);

yline(1, '-k', 'LineWidth', 1.5);
yline(-1, '-k', 'LineWidth', 1.5);
xlim([0, 5]);
ylim([-1, 1.1]);
lgd3 = legend([h_x, h_y, rx, ry], {'${\dot p_x}^\star$', '${\dot p_y}^\star$', '$\mathcal{R}_{{\dot p_x}^\star}$', '$\mathcal{R}_{{\dot p_y}^\star}$'}, 'Interpreter', 'latex', 'Location', 'best','fontsize',legendfontsize,'Position',[0.5982 0.5970 0.0763 0.2338]);
set(lgd3, 'Box', 'off', 'Color', 'none');
xlabel('time [-]', 'FontSize', labelFontSize);
set(gca,'YTickLabel', [], 'FontSize', tickFontSize);

subplot(1,4,4);
hold on
grid on;

upper_y = v_sol(1,:)+tubes_u(1,:);
lower_y = v_sol(1,:)-tubes_u(1,:);
for i = 1:length(upper_y)-1
    rectangle('Position', [(i-1)*params.dt, lower_y(i), params.dt, upper_y(i)-lower_y(i)], 'FaceColor', [color3, 0.5], 'EdgeColor', 'none');
end

upper_y = v_sol(2,:)+tubes_u(2,:);
lower_y = v_sol(2,:)-tubes_u(2,:);
for i = 1:length(upper_y)-1
    rectangle('Position', [(i-1)*params.dt, lower_y(i), params.dt, upper_y(i)-lower_y(i)], 'FaceColor', [color4, 0.5], 'EdgeColor', 'none');
end

upper_y = v_sol(1,:)+tubes_u(1,:);
lower_y = v_sol(1,:)-tubes_u(1,:);

upper_y = v_sol(2,:)+tubes_u(2,:);
lower_y = v_sol(2,:)-tubes_u(2,:);

stairs(0:params.dt:params.T, [u_cl(2,1:end-1), u_cl(2,end-1)],'k--','LineWidth', 1.1, 'MarkerSize', 10);
stairs(0:params.dt:params.T, [u_cl(1,1:end-1), u_cl(1,end-1)],'k--','LineWidth', 1.1, 'MarkerSize', 10);
stairs(0:params.dt:params.T, [u_cl_init(2,1:end-1), u_cl_init(2,end-1)],'k:', 'LineWidth', 1.1, 'MarkerSize', 10);
stairs(0:params.dt:params.T, [u_cl_init(1,1:end-1), u_cl_init(1,end-1)],'k:','LineWidth', 1.1, 'MarkerSize', 10);

h_x = stairs(0:params.dt:params.T, [v_sol(1,1:end-1), v_sol(1,end-1)],'color',color3, 'LineWidth', 1.5, 'MarkerSize', 10);
h_y = stairs(0:params.dt:params.T, [v_sol(2,1:end-1), v_sol(2,end-1)],'color',color4, 'LineWidth', 1.5, 'MarkerSize', 10);

yline(1, '-k', 'LineWidth', 1.5);
yline(1, '-k', 'LineWidth', 1.5);
xlim([0, 5]);
ylim([-1, 1.1]);
rx = fill(nan, nan,[0.8, 0.8, 0.8], 'FaceColor',color3, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
ry = fill(nan, nan,[0.8, 0.8, 0.8], 'FaceColor',color4, 'EdgeColor', 'none', 'FaceAlpha', 0.5);


lgd4 = legend([h_x, h_y, rx, ry], {'$v_x^\star$', '$v_y^\star$','$\mathcal{R}_{v_x^\star}$','$\mathcal{R}_{v_y^\star}$'}, 'Interpreter', 'latex', 'Location', 'southeast','fontsize',legendfontsize, 'Position', [0.8313 0.1049 0.0762 0.2338]);
set(lgd4, 'Box', 'off', 'Color', 'none');
xlabel('time [-]', 'FontSize', labelFontSize);

ylabel('inputs', 'FontSize', labelFontSize);
set(gca, 'FontSize', tickFontSize);

ax = axes('Position', [0 0 1 1], 'Visible', 'off');
dummy_plot = nan(1, 3);
hold on;

dummy_plot(1) = plot(ax, nan,'k:','LineWidth', 1.1, 'MarkerSize', 10);
dummy_plot(2) = plot(ax, nan,'k--','LineWidth', 1.1, 'MarkerSize', 10);
dummy_plot(3) = plot(ax, nan, '-k', 'LineWidth', 1.5);

color_labels = { 'non-robust', 'sample','$\mathcal{C}$'};
hold off;
lgd = legend(ax, color_labels, 'Location', 'southeastoutside','Interpreter','latex','fontsize',legendfontsize,'Position',[0.9068 0.1746 0.1061 0.1402]);

set(lgd,'Box', 'off', 'Color', 'none');

width_cm = 25; % Width in centimeters
height_cm = 8; % Height in centimeters
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 width_cm height_cm]);

% filename = './figs/fig1';
% exportgraphics(gcf, [filename,'.pdf'], 'ContentType', 'vector', 'Resolution', 1200);
% saveas(gcf, [filename,'.png']);

%% Figure 2
figure(2);

custom_colors = [0.0504    0.0298    0.5280
                 0.6107    0.0902    0.6200
                 0.9283    0.4730    0.3261
                 0.9400    0.9752    0.1313];

custom_colors(end,:) = [];
subplot_titles = {'$p_x^\star$', '$p_y^\star$', '$\psi^\star$', '${\dot p_x}^\star$', '${\dot p_y}^\star$', '$\dot \psi^\star$'};

labelFontSize = 10;
tickFontSize = 10;
B_stack = cat(3,B2,B3-B2-B1,B1);

for i = 1:6
    s = subplot(1,6,i);
    b = bar(squeeze(B_stack(i,:,:)), 'stacked');
    hold on;
    for j = 1:size(custom_colors,1)
        b(j).FaceColor = custom_colors(j,:);
    end
    grid on;
    ylim([0.0005, 0.1]);
    
    xlabel('step $k$', 'FontSize', labelFontSize,'Interpreter','latex');
    title(subplot_titles{i},'Interpreter','latex', 'FontSize', labelFontSize);
end
color_labels = {'$\mathcal{W}$', '${\Theta}^\star$', '$\tau^{\star 2} \mu$'};
subplot(1,6,1);
ylabel('$\sigma_{i,k}^\star$', 'FontSize', labelFontSize,'Interpreter','latex');
    set(gca, 'YScale', 'log', 'XTickLabel', [],'FontSize', tickFontSize);
    
for i=2:6
    subplot(1,6,i);    
    set(gca, 'YTickLabel', [], 'XTickLabel', [],'YScale', 'log','FontSize', tickFontSize);
end
    
ax = axes('Position', [0 0 1 1], 'Visible', 'off');
dummy_bars = nan(1, numel(color_labels));
hold on;
for i = 1:numel(color_labels)
    dummy_bars(i) = bar(ax, nan, 'FaceColor', custom_colors(i, :));
end
hold off;
lgd = legend(ax, color_labels, 'Location', 'eastoutside','Interpreter','latex');
set(lgd, 'Box', 'off', 'Color', 'none');
ax.XAxis.Color = 'none';
ax.YAxis.Color = 'none';

width_cm = 25; % Width in centimeters
height_cm = 5; % Height in centimeters
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 width_cm height_cm]);

% filename = './figs/fig2';
% exportgraphics(gcf, [filename,'.pdf'], 'ContentType', 'vector', 'Resolution', 1200);
% saveas(gcf, [filename,'.png']);