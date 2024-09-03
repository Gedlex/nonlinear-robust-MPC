function [solver,bx,bg] = setupMPC(obj,terminal_constraint,solver_options)
    arguments
        obj (1,1) CCM 
        terminal_constraint (1,1) logical = true
        solver_options (1,1) struct = struct('ipopt', struct('max_iter',500, 'print_level',0), 'print_time',0)
    end
    import casadi.*
    
    % Compute number of variables / parameters
    n_var = obj.n_v + obj.n_z;      % number of optimization vars
    n_par = obj.nx;                 % number of parameters
    
    % Create optimization variables / parameters
    y = MX.sym('y', n_var);
    x_t = MX.sym('par',n_par);
    
    % Extract variables
    v = reshape(y(1:obj.n_v),             [obj.nu,obj.N]);
    z = reshape(y(obj.n_v + (1:obj.n_z)), [obj.nx,obj.N+1]);
    
    % Get reference state and input that keep the nominal system in equilibrium
    z_ref = obj.params.x_ref;
    v_ref = obj.params.u_ref(z_ref,zeros(obj.np,1));
    
    %% Compute Cost
    objective = getObjective(obj,z,v,z_ref,v_ref);
    
    %% Equality constraints
    ceq=[];
    cineq=[];
    % Dynamic constraints
    for i = 1:obj.N
        ceq = [ceq; z(:,i+1) - obj.sys.ddyn(z(:,i),v(:,i))];
    end
    % Constraint on initial state
    ceq = [ceq; z(:,1) - x_t];
    
    % Terminal equality constraint
    if terminal_constraint
        ceq = [ceq; z(:,end) - z_ref];
    end
    
    %% Inequality constraints
    for k = 1:obj.N
        % State constraints
        cineq = [cineq; obj.params.F_x*z(:,k) - obj.params.b_x];
        
        % Input constraints
        cineq = [cineq; obj.params.F_u*v(:,k) - obj.params.b_u];
        
        % Obstacle constraints
        if isa(obj.params.h_obs,"function_handle")
            cineq = [cineq; obj.params.h_obs(z(:,k))];
        end
    end
    % Terminal state constraints
    cineq = [cineq; obj.params.F_x*z(:,end) - obj.params.b_x];
    
    %% Setup NPSOL
    n_eq = length(ceq);
    n_ineq = length(cineq);
    
    % Define lower and upper bound on constraints
    bg.clb = [zeros(n_eq,1); -inf(n_ineq,1)];
    bg.cub =  zeros(n_eq     +    n_ineq,1);
    
    % Define lower and upper bound on geodesic variables
    bx.ylb =-inf(n_var,1);
    bx.yub = inf(n_var,1);
    
    % Initialize optimizer
    con = [ceq;cineq];
    nlp = struct('x', y, 'p', x_t, 'f', objective, 'g', con);
    solver = nlpsol('solver', 'ipopt', nlp, solver_options);
end