function params = CSparams1()
% CSPARAMS1 initializes a struct containing the controller parameters for 
% a post capture stabilization.
    
    params = Parameters();
    params.label = 'Capture Stabilization';
    params.version = 1;
    
    % Horizon and discrete step time
    params.N = 10;
    params.dt = 0.5;
    
    % Reference state and input
    params.x_ref = zeros(6,1);
    params.u_ref = @(x,theta) [0;0];
    params.x0 = [.7;.7;.5;.5;.5;.5];
    
    % Constraints and vertices
    params.F_u = [eye(2); -eye(2)];
    params.b_u = ones(2*2,1);
    params.F_x = [eye(6); -eye(6)];
    params.b_x = ones(2*6,1);
    params.theta_v = [-0.01, 0.01];
    params.theta_true = 0.01;
    params.w_max = 0.002;
    params.nw = 3;
    
    % Cost
    params.Q_cost = eye(6);
    params.R_cost = eye(2);
    params.P_cost = 30*params.Q_cost;

    % Integrator
    params.integrator = "multi";
end