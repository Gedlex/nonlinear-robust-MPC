%% Test Cart-Pole Dynamics
clearvars;
close all;
clc;

%% Solve SLS
% Load parameters
params = CPparams1();

% Create pre-stabilized dynamics
sys = CartPole(params.dt,integrator=params.integrator,param_uncertainty=params.param_uncertainty);

% Create sls class
sls = SLS(sys,params);

%% Solve Problem
[v, s, y0] = sls.solve(params.x0);

% Notify user
Data = load('gong.mat');
sound(Data.y, Data.Fs)

%% Export and Plot Data
% Define axis limits
limits_x = {[-5,       5];
            [-6*pi, 6*pi];
            [-5,       5];
            [-6*pi, 6*pi]};
limits_u = {[];
            []};

% Generate output path
out_path = sprintf('%s_N=%d_dt=%.3f_wmax=%.3f_theta=[%.4f_%.4f]_int=%s_max_iter=%d', ...
                   class(sys),params.N,params.dt,params.w_max,params.theta_v(1),params.theta_v(2),sys.integrator,10000);

% Save data
save(['./data/SLS_',out_path,'.mat'],'-struct','s')

% Plot trajectories
sys.plot.trajectories(s.t,s.z_sol,s.v_sol(:,1:params.N),s.tubes_x,s.tubes_u(:,1:params.N),'SLS',xLimits=limits_x,uLimits=limits_u,sharedAxes=false);
% sys.plot.exportPlot(gcf, ['./figs/SLS_',out_path,'.pdf'], 'pdf')