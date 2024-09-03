%% Compare results
clc, clearvars, close all

%% Load data
% sls_data = importdata('data/SLS_PlanarQuadrotor_N=10_dt=0.075_wmax=0.2000_theta=[-0.00_0.00]_int=rk4.mat');
% ccm_data = importdata('data/CCM_PlanarQuadrotor_N=10_dt=0.075_wmax=0.2000_theta=[0.00_0.00]_int=rk4.mat');
ccm_data = importdata('data/CS/CCM_CaptureStabilization_N=10_dt=0.500_wmax=0.010_theta=[0.0000_0.0000]_int=rk4_max_iter=10000_rho=-0.200.mat');

%% Create Plots
plt = PlanarQuadrotorPlots();

% Define axis limits
limits_x = {[-1, 1];
            [-1, 1];
            [-1, 1];
            [-1, 1];
            [-1, 1];
            [-1, 1]};
limits_u = {[-1, 1];
            [-1, 1]};

% Generate output path
out_path = sprintf('Planar_Quadrotor_s_N=%d_dt=%.3f_wmax=%.4f_int=%s',10,0.075,0.2,'rk4');

% Plot trajectories
% plt.trajectories(sls_data.t,sls_data.z_sol,sls_data.v_sol(:,1:end-1)sls_data.tubes_x,sls_data.tubes_u(:,1:end-1),'SLS', ...
%                  xLimits=limits_x,uLimits=limits_u,alphaArea=0.3);
plt.trajectories(ccm_data.t,ccm_data.z_sol,ccm_data.v_sol,ccm_data.tubes_x,ccm_data.tubes_u,'CCM', ...
                 xLimits=limits_x,uLimits=limits_u,alphaArea=0.3);
%plt.exportPlot(gcf, ['./figs/CCM_SLS_traj_',out_path,'.svg'])