%% Main
% -------------------------------------------------------------------------
% File: main.m
% Author: Alexander Erdin (aerdin@ethz.ch)
% Date: 22th May 2024
% License: MIT
% Reference:
%
% -------------------------------------------------------------------------
clearvars;
close all;
clc;

%% Compute Prestabilizing Feedback
% Load parameters
use_sos=true;
params = PQRparams3();

% Create planar quadrotor
sys = PlanarQuadrotor(params.dt,integrator=params.integrator,param_uncertainty=params.param_uncertainty);

% Create ccm class instance with pre-stabilizing feedback
ccm = CCM(sys,params, ...
          use_sos=use_sos, ...
          monomial_degree=0, ...
          terminal_constraint=false);

% Get pre-stabilizing feedback matrix K
K = ccm.Y_fcn()/ccm.W_fcn();

% Solve problem
% [v,s,y0] = ccm.solve(params.x0);

%% Initialize Prestabilized SLS
% Load parameters
params = PQRparams2();
params.integrator = "rkdp";
params.dt = 0.1;

% Create pre-stabilized dynamics
sys = PlanarQuadrotorPreStabilized(params.dt,K=K,integrator=params.integrator,param_uncertainty=params.param_uncertainty);

% Create sls class
sls = SLS(sys,params,K=K)%recompute_mu=true);

%% Initialize Regular SLS
% Load parameters
params = PQRparams2();
params.integrator = "rkdp";
params.dt = 0.1;

% Create pre-stabilized dynamics
sys = PlanarQuadrotor(params.dt,integrator=params.integrator,param_uncertainty=params.param_uncertainty);

% Create sls class
sls = SLS(sys,params);

%% Solve Problem
[v, s, y0] = sls.solve(params.x0);

% Notify user
Data = load('gong.mat');
sound(Data.y, Data.Fs)

%% Export and Plot Data
% Define axis limits
limits_x = {[];
            [];
            [-pi/3, pi/3];
            [-2,       2];
            [-1,       1];
            [-pi/1,     pi/1]};
limits_u = {[-1,     3.5];
            [-1,     3.5]};
% limits_x = repmat({[-1, 1]},sys.nx,1);
% limits_u = repmat({[-1, 1]},sys.nu,1);

% Generate output path
out_path = sprintf('%s_N=%d_dt=%.3f_wmax=%.3f_theta=[%.4f_%.4f]_int=%s_max_iter=%d', ...
                   class(sys),params.N,params.dt,params.w_max,params.theta_v(1),params.theta_v(2),sys.integrator,10000);

% Save data
save(['./data/SLS_',out_path,'.mat'],'-struct','s')

% Plot trajectories
sys.plot.trajectories(s.t,s.z_sol,s.v_sol(:,1:params.N),s.tubes_x,s.tubes_u(:,1:params.N),'SLS',xLimits=limits_x,uLimits=limits_u,sharedAxes=false);
sys.plot.exportPlot(gcf, ['./figs/SLS_',out_path,'.pdf'], 'pdf')

% Plot errors
% sys.plot.errors(s.B2,s.B3-s.B2-s.B1,s.B1)
% sys.plot.exportPlot(gcf, ['./figs/errors_',out_path,'.jpeg'], 'jpeg')

% Plot xy trajectory
% sys.plot.xy(s.z_sol(1,:),s.z_sol(2,:),'SLS',limits={[-2,2];[-2,2]})
% sys.plot.exportPlot(gcf, ['./figs/xy_',out_path,'.jpeg'], 'jpeg')

%% Plots from original paper
if strcmp('Yes', questdlg('Continue to plot figures from original paper?', 'Plot Figures', 'Yes','No','No'))
    load(['./data/SLS_',out_path,'.mat'])
    plot_paper;
end