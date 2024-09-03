clearvars;
clc;

% Load parameters
use_sos = false;
params = CPparams1();

% Create planar quadrotor
sys = CaptureStabilization(params.dt, ...
                           integrator=params.integrator, ...
                           approximate=use_sos, ...
                           param_uncertainty=params.param_uncertainty);

% Define rhos
rho = -0.1:0.01:-0.01;

fprintf('Job started.\n')
for i = 1:numel(rho)
    fprintf('Solve CCM for rho = %.3f\n',rho(i))
    try
        % Create ccm class instance
        ccm = CCM(sys,params, ...
                  rho=rho(i), ...
                  W_idxs=3, ...
                  use_sos=use_sos, ...
                  monomial_degree=14, ...
                  terminal_constraint=false);
        
        % Solve problem
        [v, s, y0] = ccm.solve(params.x0);
        
    catch ME
        disp(getReport(ME,'extended'));
        continue
    end
    
    % Generate output path
    out_path = sprintf('%s_N=%d_dt=%.3f_wmax=%.3f_theta=[%.4f_%.4f]_int=%s_max_iter=%d_rho=%.3f', ...
                        class(sys),params.N,params.dt,params.w_max,params.theta_v(1),params.theta_v(2),sys.integrator,10000,ccm.rho_des);
    
    % Save data
    save(['./data/CCM_',out_path,'.mat'],'-struct','s')
end

fprintf('Job finished!\n')