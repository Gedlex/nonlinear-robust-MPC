clearvars;
clc;

% Load parameters
params = PQRparams3();

% Define rhos
rho = 0.1:0.1:0.9;

% Create planar quadrotor
sys = PlanarQuadrotor(params.dt,integrator=params.integrator,param_uncertainty=params.param_uncertainty);

% Create ccm class instance with pre-stabilizing feedback
ccm = CCM(sys,params, ...
          rho=rho, ...
          use_sos=true, ...
          monomial_degree=0, ...
          terminal_constraint=false);

% Get pre-stabilizing feedback matrix K
K = ccm.Y_fcn()/ccm.W_fcn();

% Create pre-stabilized planar quadrotor
sys = PlanarQuadrotorPreStabilized(params.dt,K=K,integrator=params.integrator,param_uncertainty=params.param_uncertainty);

% Define horizons
N = 1:30;

fprintf('Job started.\n')
for i = 1:numel(N)
    % Set horizon
    params.N = N(i);

    fprintf('Solve SLS for N = %d\n',N(i))
    try
        % Create sls class instance
        sls = SLS(sys,params,K=K);

        % Solve problem
        [v, s, y0] = sls.solve(params.x0);

    catch ME
        disp(getReport(ME,'extended'));
        continue
    end

    % Generate output path
    out_path = sprintf('%s_N=%d_dt=%.3f_wmax=%.3f_theta=[%.4f_%.4f]_int=%s_max_iter=%d', ...
                       class(sys),params.N,params.dt,params.w_max,params.theta_v(1),params.theta_v(2),sys.integrator,10000);

    % Save data
    save(['./data/PQR2/SLSK_',out_path,'.mat'],'-struct','s')
end

fprintf('Job finished!\n')