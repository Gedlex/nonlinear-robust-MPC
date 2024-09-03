clearvars;
clc;

% Load parameters
params = PQRparams2();

% Create planar quadrotor
sys = PlanarQuadrotor(params.dt,integrator=params.integrator,param_uncertainty=params.param_uncertainty);

% Define disturbances
w_max = linspace(0.1,0.3,11);

fprintf('Job started.\n')
for i = 1:numel(w_max)
    % Set disturbance
    params.w_max = w_max(i);

    fprintf('Solve SLS for w_max = %.3f\n',w_max(i))
    try
        % Create sls class instance
        sls = SLS(sys,params);

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
    save(['./data/PQR/SLS_',out_path,'.mat'],'-struct','s')
end

fprintf('Job finished!\n')