function params = PQRparams3()
% PQRPARAMS3 initializes a struct containing the controller parameters for 
% a planar quadrotor.
    
    params = Parameters();
    params.label = 'Planar Quadrotor';
    params.version = 1;
    
    % Horizon and discrete step time
    params.N = 10;
    params.dt = 0.075;
    
    % Reference state and input
    params.x_ref = zeros(6,1);
    params.u_ref = @(x,theta) zeros(2,1);
    params.x0 = [-2;-15;0;0;0;0];
    
    % Constraints and vertices
    params.F_u = [eye(2); -eye(2)];
    params.b_u = [3.5,3.5,1,1]';
    params.F_x = blkdiag([eye(2); -eye(2)], [eye(4); -eye(4)]);
    params.b_x = [[20, 20, 20, 20],         [pi/20, 2, 1, pi/10, pi/20, 2, 1, pi/10]]';
    params.theta_v = zeros(1,2);
    params.theta_true = 0;
    params.w_max = 0.1;
    params.nw = 1;
    
    % Cost
    params.Q_cost = eye(6);
    params.R_cost = eye(2);
    params.P_cost = 30*params.Q_cost;
    params.regularizer = 1E-6;
    
    % Integrator
    params.integrator = "rkdp";
end