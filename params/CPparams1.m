function params = CPparams1()
% CPPARAMS1 initializes a struct containing the controller parameters for 
% a cart-pole system.
    
    params = Parameters();
    params.label = 'Cart-Pole';
    params.version = 1;
    
    % Horizon and discrete step time
    params.N = 10;
    params.dt = 0.075;
    
    % Reference state and input
    params.x_ref = zeros(4,1);
    params.u_ref = @(x,theta) 0;
    params.x0 = [-1;pi/6;0;0];
    
    % Constraints and vertices
    params.F_u = [1; -1];
    params.b_u = [100, 100]';
    params.F_x = [eye(4); -eye(4)];
    params.b_x = [5, 6*pi, 5, 6*pi, 5, 6*pi, 5, 6*pi]';
    params.theta_v = zeros(1,2);
    params.theta_true = 0;
    params.w_max = 0.05;
    params.nw = 2;
    
    % Cost
    params.Q_cost = eye(4);
    params.R_cost = 0.01*eye(1);
    params.P_cost = [1,  0,  0,  0;
                     0,  30, 0,  0;
                     0,  0,  1,  0;
                     0,  0,  0,  1];
    params.regularizer = 1E-6;
    
    % Integrator
    params.integrator = "rk4";
end