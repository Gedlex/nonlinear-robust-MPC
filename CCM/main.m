%% Main
% -------------------------------------------------------------------------
% File: main.m
% Author: Alexander Erdin (aerdin@ethz.ch)
% Date: 19th April 2024
% License: MIT
% Reference:
%
% -------------------------------------------------------------------------
clearvars;
close all;
clc;

% Load parameters
use_sos = false;
params = CPparams3();

% Create planar quadrotor
sys = CaptureStabilization(params.dt, ...
                           integrator=params.integrator, ...
                           approximate=use_sos, ...
                           param_uncertainty=params.param_uncertainty);

% Create ccm class instance
ccm = CCM(sys,params, ...
          rho=-0.3, ...
          W_idxs=3, ...
          use_sos=use_sos, ...
          monomial_degree=14, ...
          terminal_constraint=false);

%% Solve problem
[v,s,y0] = ccm.solve(params.x0);

%% Export and Plot Data
% Define axis limits
% limits_x = {[];
%             [];
%             [-pi/3, pi/3];
%             [-2,       2];
%             [-1,       1];
%             [-pi,     pi]};
% limits_u = {[-1,     3.5];
%             [-1,     3.5]};
limits_x = {[-1, 1];
            [-1, 1];
            [-1, 1];
            [-1, 1];
            [-1, 1];
            [-1, 1]};
limits_u = {[-1, 1];
            [-1, 1]};

% Generate output path
out_path = sprintf('%s_N=%d_dt=%.3f_wmax=%.4f_theta=[%.2f_%.2f]_int=%s', class(sys),params.N,params.dt,params.w_max,params.theta_v(1),params.theta_v(2),sys.integrator);

% Save data
s.limits_x = limits_x;
s.limits_u = limits_u;
% save(['./data/CCM_',out_path,'.mat'],'-struct','s')

% Plot trajectories
sys.plot.trajectories(s.t,s.z_sol,s.v_sol,s.tubes_x,s.tubes_u,'CCM',xLimits=limits_x,uLimits=limits_u);

%% Closed Loop Simulation (From Paper)
% Number of simulation steps
N_sim = 30;

% Load parameters
adaptive = false;
exportname = ''; % 'res_robust_half'; %'res_adaptive'; %'res_robust'; %
[params, obs] = PQRparams_theta(); %PQRparams_theta_half(); %

% Create planar quadrotor
sys = PlanarQuadrotor(params.dt, ...
                      integrator=params.integrator, ...
                      param_uncertainty=params.param_uncertainty);

% Create ccm class instance
ccm = CCM(sys,params, ...
          adaptive=adaptive, ...
          exact_initialization=false, ...
          terminal_constraint=false, ...
          use_sos=true);

%%

% Compute initial solver guess (for non-adpative RMPC)
if ~adaptive && ~contains(exportname, 'half')
    % Initial guess for states and inputs (first to the midpoint, then to origin)
    mp = [-0.5;-2; zeros(4,1)];
    N_mp = round((ccm.N+1)/3);
    z0 = zeros(ccm.nx,ccm.N+1);
    z0(:,1:N_mp)   = ccm.params.x0 + (mp - ccm.params.x0)*linspace(0,1,N_mp);
    z0(:,N_mp:end) = mp + (ccm.params.x_ref - mp)*linspace(0,1,ccm.N+2-N_mp);
    v0 = ccm.params.u_ref(ccm.params.x_ref,zeros(ccm.np,1));
    
    idx = ccm.n_geod;
    y0 = zeros(ccm.n_var,1);
    y0(idx + (1:ccm.n_v)) = repmat(v0,ccm.N,1);             idx = idx+ccm.n_v;
    y0(idx + (1:ccm.n_z)) = reshape(z0,ccm.nx*(ccm.N+1),1); idx = idx+ccm.n_z;
    
    % Set initial solver value
    ccm.y0 = y0;
end

% Initialize parameters
x = ccm.params.x0;
Theta = ccm.params.theta_v;

% Solve closed loop RAMPC
tic;
sol_store = cell(N_sim,1);
for k = 1:N_sim
    fprintf('Simulation step %d/%d\n',k,N_sim);
    % Solve mpc
    [v_sol,s,y0] = ccm.solve(x,Theta);

    % Update initial solver value
    ccm.y0 = y0;

    % Simulate real dynamics using rk4 (rk4 steps [t_k, t_k+dt/2, t_k+dt])
    z_rk = zeros(ccm.nx,3);
    z_rk(:,1) = s.z_sol(:,1);
    z_rk(:,3) = s.z_sol(:,2);
	z_rk(:,2) = ccm.sys.ddyn(s.z_sol(:,1),v_sol,theta=s.theta_bar_sol,dt=ccm.params.dt/2);
    v_rk = repmat(v_sol,1,3);                  % Input is piece-wise constant
    d_rk = ccm.params.w_max*(2*rand(1,3) - 1); % Generate random disturbances
    
    % Compute feedback kappa based on nominal state and input
    f = @(x,u,d) ccm.sys.fw(x,u,ccm.params.theta_true,d);
    kappa_rk = @(x,z,v) ccm.kappa(x,z,v);
    [x_new, u] = dynamics_real_RK(f,x,kappa_rk,z_rk,v_rk,d_rk,ccm.params.dt);
    
    % Perform set membership update
    if adaptive 
        % Get noisy measurements for set-membership updates
        dx_noise = ccm.params.w_max*(2*rand(ccm.nx,3) - 1); %#ok
        meas_dx = ccm.sys.fw(x,u,ccm.params.theta_true,d_rk(1)) + dx_noise;
        
        % Set membership estimation
        non_fals = compute_non_falsified(x,u,meas_dx,ccm.params.w_max,ccm.params.w_max,ccm.sys.m_nom);
        Theta = intersection(Theta,non_fals);
    end
    
    % Store solution
    s.x = x;
    s.u = u;
    sol_store{k} = s;
    
    % Update x
    x = x_new;
end

% Get results
sol_store = cell2mat(sol_store);
res.z = cat(3,sol_store.z_sol);
res.v = cat(3,sol_store.v_sol);
res.traj = [sol_store.x];
res.delta_tight = [sol_store.delta_tight];
res.M_pos = ccm.M_under(1:2,1:2);
res.N = ccm.N;
res.obs = obs;

% Plot results
figure(); hold on;
plot(res.traj(1,:), res.traj(2,:));
plot(res.z(1,:,1), res.z(2,:,1),'+-');

% Plot tube
for j = 1:res.N
    plot_tube(res.z(:,j,1),res.M_pos,res.delta_tight(j,1),[0,0,1,0.5],1);
end

% Plot obstacles
visualize_obs(obs);

% Save results
if ~isempty(exportname); save(['data/sim_results/', exportname, '.mat'],'-struct','res'); end