function [params, obs] = PQRparams_theta()
% PQRPARAMS_THETA initializes a struct containing the controller parameters
% for a planar quadrotor.
    
    params = Parameters();
    params.label = 'Planar Quadrotor';
    params.version = 1;
    
    % Horizon and discrete step time
    params.N = 25;
    params.dt = 0.15;
    
    % Reference state and input
    params.x_ref = zeros(6,1);
    params.u_ref = @(x,theta) ones(2,1)*9.81/(2*(2.0576 + theta(:)));
    params.x0 = [-2;-2;0;0;0;0];
    
    % Constraints and vertices
    params.F_u = [eye(2); -eye(2)];
    params.b_u = [3.5,3.5,1,1]';
    params.F_x = blkdiag([eye(2); -eye(2)], [eye(4); -eye(4)]);
    params.b_x = [[20, 20, 20, 20],         [pi/3, 2, 1, pi, pi/3, 2, 1, pi]]';
    params.theta_v = 2.0576*[-0.01, 0.01];
    params.theta_true = 0.0206;
    params.w_max = 0.1;
    params.nw = 1;

    % Obstacle constraints
    obs = [-1.5   -0.9    0.16;
           -1.0   -0.5    0.16;
           -0.7   -1.2    0.16;
            0.3   -1.0    0.16];
    params.h_obs = @(x) -sqrt((x(1) - obs(:,1)).^2 + (x(2) - obs(:,2)).^2) + obs(:,3);
    
    % Cost
    params.Q_cost = eye(6);
    params.R_cost = 0.1*eye(2);
    params.P_cost = zeros(6);
    params.regularizer = 1E-6;
    
    % Integrator
    params.integrator = "rk4";
end