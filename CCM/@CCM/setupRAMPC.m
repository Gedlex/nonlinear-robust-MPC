function obj = setupRAMPC(obj,terminal_constraint,exact_initialization,solver_options)
    arguments
        obj (1,1) CCM 
        terminal_constraint (1,1) logical = true
        exact_initialization (1,1) logical = false
        solver_options (1,1) struct = struct('ipopt', struct('max_iter',500, 'print_level',0), 'print_time',0)
    end
    import casadi.*
    
    % Compute number of variables / parameters
    n_var = obj.n_geod + obj.n_v + obj.n_z + obj.n_delta + obj.n_w_bar; % number of optimization vars
    n_par = obj.nx + 2*obj.np;                                          % number of parameters
    if obj.adaptive
        n_var = n_var + obj.np + 1;
    end
    
    % Create optimization variables / parameters
    y = MX.sym('y', n_var);
    par = MX.sym('par',n_par);
    
    % Extract variables
    idx = 0;
    c_geod = y(1:obj.n_geod);                             idx = idx+obj.n_geod;
    v = reshape(y(idx+(1:obj.n_v)), [obj.nu,obj.N]);      idx = idx+obj.n_v;
    z = reshape(y(idx+(1:obj.n_z)), [obj.nx,obj.N+1]);    idx = idx+obj.n_z;
    delta = y(idx+(1:obj.n_delta));                       idx = idx+obj.n_delta;
    w_bar = reshape(y(idx+(1:obj.n_w_bar)), [4,obj.N]);   idx = idx+obj.n_w_bar;
    if obj.adaptive
        theta_bar = y(idx+(1:obj.np));                    idx = idx+obj.np;
        delta_bar = y(idx+1);
    else
        theta_bar = zeros(obj.np,1);
        delta_bar = min(obj.delta_bar_x,obj.delta_bar_u);
    end
    
    % Extract parameters
    x_t = par(1:obj.nx);
    Theta_t = reshape(par(obj.nx+1:end),2,obj.np)';  % [min, max]
    
    % Get reference state and input that keep the nominal system in equilibrium
    z_ref = obj.params.x_ref;
    v_ref = obj.params.u_ref(z_ref,theta_bar);
    
    % Get all vertex combinations of d and theta
    Theta_vert = getCombinations(repmat([1,2],obj.np,1)); % workaround to deal with MX variable
    d_vert = getCombinations(obj.params.disturbance_v);
    
    % Create function for CCM
    xc = SX.sym('x',obj.nx);
    Wc = obj.W_fcn(xc);
    W_chol = Function('W_chol',{xc},{chol(Wc)});
    
    %% Compute Cost
    objective = getObjective(obj,z,v,z_ref,v_ref);
    
    % Add regularizer
    objective = objective + obj.params.regularizer*(delta'*delta);
    
    %% Compute Discretized Dynamics
    % Initialize constraints
    ceq=[];
    cineq=[]; % all inequalty constraints are formulated as cineq <= 0
    for i = 1:obj.N        
        % Dynamics - Runge-Kutta 4 discretization; input is piece-wise constant could instead also treat k1,...,k4 as decision variables
        k1=obj.sys.fw(z(:,i),v(:,i),theta_bar,0);
        k2=obj.sys.fw(z(:,i)+obj.sys.dt/2*k1,v(:,i),theta_bar,0);
        k3=obj.sys.fw(z(:,i)+obj.sys.dt/2*k2,v(:,i),theta_bar,0);
        k4=obj.sys.fw(z(:,i)+obj.sys.dt*k3,v(:,i),theta_bar,0);
        ceq = [ceq; z(:,i+1) - z(:,i) - obj.sys.dt/6*(k1+2*k2+2*k3+k4)];
        
        % Tube dynamics - Runge-Kutta 4 discretization using w_bar as upper bound for max term
        k1_tube=f_delta_cont(delta(i),obj.rho,obj.L_D,w_bar(1,i));
        k2_tube=f_delta_cont(delta(i)+obj.sys.dt/2*k1_tube,obj.rho,obj.L_D,w_bar(2,i));
        k3_tube=f_delta_cont(delta(i)+obj.sys.dt/2*k2_tube,obj.rho,obj.L_D,w_bar(3,i));
        k4_tube=f_delta_cont(delta(i)+obj.sys.dt*k3_tube,obj.rho,obj.L_D,w_bar(4,i));
        ceq = [ceq; delta(i+1) - delta(i) - obj.sys.dt/6*(k1_tube+2*k2_tube+2*k3_tube+k4_tube)];
        
        % Array of intermediate RK points (used for continous-time constraints on w_bar)
        z_RK = [z(:,i), z(:,i) + obj.sys.dt/2*k1, z(:,i) + obj.sys.dt/2*k2, z(:,i) + obj.sys.dt*k3];
        delta_RK = [delta(i), delta(i) + obj.sys.dt/2*k1_tube, delta(i) + obj.sys.dt/2*k2_tube, delta(i) + obj.sys.dt*k3_tube];
        
        % Enforce w_bar >= L_g*|theta^i - theta_bar|*delta + ||G*(theta^i - theta_bar) + E*d^j||_M(z) (Eqn. (20))
        % Loop over intermediate RK points
        for j = 1:4
            % Loop over vertex combinations of theta
            for k = 1:size(Theta_vert,2)
                idx = sub2ind([obj.np, 2], (1:obj.np)', Theta_vert(:,k));
                Theta_k = Theta_t(idx);
                
                % Loop over vertex combinations of d
                for l = 1:size(d_vert,2)
                    d_l = d_vert(:,l);
                    
                    % Use two inequalities instead of abs
                    cineq = [cineq; -w_bar(j,i) + obj.L_G*(Theta_k - theta_bar)*delta_RK(j) + ...
                             obj.w_tilde(z_RK(:,j),v(:,i),theta_bar,Theta_k,d_l,W_chol(z_RK(:,j)))];
                    cineq = [cineq; -w_bar(j,i) - obj.L_G*(Theta_k - theta_bar)*delta_RK(j) + ...
                             obj.w_tilde(z_RK(:,j),v(:,i),theta_bar,Theta_k,d_l,W_chol(z_RK(:,j)))];
                end
            end
        end
    end
    
    %% Constraints
    xineq_idxs = [];
    uineq_idxs = [];
    for k = 1:obj.N
        % State constraints
        xineq_idxs = [xineq_idxs; numel(cineq) + (1:numel(obj.params.b_x))'];
        cineq = [cineq; obj.params.F_x*z(:,k) - obj.params.b_x + obj.c_x * delta(k)];
        
        % Input constraints
        uineq_idxs = [uineq_idxs; numel(cineq) + (1:numel(obj.params.b_u))'];
        cineq = [cineq; obj.params.F_u*v(:,k) - obj.params.b_u + obj.c_u * delta(k)];
        
        % Obstacle constraints
        if isa(obj.params.h_obs,"function_handle")
            cineq = [cineq; obj.params.h_obs(z(:,k)) + obj.c_obs * delta(k)];
        end
    end
    % Terminal state constraints
    xineq_idxs = [xineq_idxs; numel(cineq) + (1:numel(obj.params.b_x))'];
    cineq = [cineq; obj.params.F_x*z(:,end) - obj.params.b_x + obj.c_x * delta(end)];
    
    %% Theta_bar is in Theta_0
    if obj.adaptive
        cineq = [cineq;  obj.params.theta_v(:,1) - theta_bar];
        cineq = [cineq; -obj.params.theta_v(:,2) + theta_bar];
    end
    
    %% Initial Tube Scaling Constraint (V_delta <= delta_0)
    % Geodesic equality constraints
    ceq = [ceq; obj.geod.Aeq*c_geod - [z(:,1); x_t]];
    
    % Delta_0 >= argmax_{gamma,gamma_s} V_delta(gamma,gamma_s,x,z)
    gamma = reshape(c_geod,obj.geod.D+1,obj.nx)'*obj.geod.T;
    gamma_s = reshape(c_geod,obj.geod.D+1,obj.nx)'*obj.geod.Tdot;
    cineq = [cineq; -delta(1)^2 + VdeltaSquared(obj.geod.N,gamma,gamma_s,obj.geod.w_cheby,obj.W_fcn)];
    
    % Ensure that delta_0 == 0 (for comparison with other methods)
    if exact_initialization
        ceq = [ceq; delta(1)];
    else % delta_0 >= 0 (delta_i >= 0 follows for i > 0)
        cineq = [cineq; -delta(1)];
    end
    
    % Ensure that gamma lies in state constraints for faster convergence (theoretically, this is already guaranteed by Proposition 5)
    for i = 1:size(gamma,2)
        cineq = [cineq; obj.params.F_x*gamma(:,i) - obj.params.b_x];
    end
    
    %% Terminal equality constraint (steady state does not depend on theta_bar)
    if terminal_constraint
        ceq = [ceq; z(:,end) - z_ref];
        cineq = [cineq; delta(end) - delta_bar];
    end
    
    if obj.adaptive
        % Delta_bar is largest value for which the thightened constraints are satisfied in the terminal set
        cineq = [cineq; delta_bar - obj.delta_bar_x]; % steady-state does not depend on theta_bar
        
        cineq = [cineq; obj.params.F_u*v_ref - obj.params.b_u + obj.c_u * delta_bar];
        
        % Tube dynamics must be invariant for delta_bar
        % Loop over vertex combinations of theta
        for k = 1:size(Theta_vert,2)
            idx = sub2ind([obj.np, 2], (1:obj.np)', Theta_vert(:,k));
            Theta_k = Theta_t(idx);
            
            % Loop over vertex combinations of d
            for l = 1:size(d_vert,2)
                d_l = d_vert(:,l);
                
                % Use two inequalities instead of abs
                cineq = [cineq; -(obj.rho - obj.L_D)*delta_bar + ...
                         obj.L_G*(Theta_k - theta_bar)*delta_bar + ...
                         obj.w_tilde(z_ref,v_ref,theta_bar,Theta_k,d_l,W_chol(z_ref))];
                cineq = [cineq; -(obj.rho - obj.L_D)*delta_bar + ...
                         obj.L_G*(-Theta_k + theta_bar)*delta_bar + ...
                         obj.w_tilde(z_ref,v_ref,theta_bar,Theta_k,d_l,W_chol(z_ref))];
            end
        end
    end
    
    %% Setup NPSOL
    obj.n_eq = length(ceq);
    obj.n_ineq = length(cineq);
    obj.n_var = n_var;
    
    % Shift idxs by equality constraints
    obj.xineq_idxs = xineq_idxs + numel(ceq);
    obj.uineq_idxs = uineq_idxs + numel(ceq);
    
    % Define lower and upper bound on constraints
    obj.clb = [zeros(obj.n_eq,1); -inf(obj.n_ineq,1)];
    obj.cub =  zeros(obj.n_eq     +    obj.n_ineq,1);
    
    % Define lower and upper bound on geodesic variables
    obj.ylb =-inf(n_var,1);
    obj.yub = inf(n_var,1);
    obj.ylb(1:obj.n_geod) = obj.geod.lb;
    obj.yub(1:obj.n_geod) = obj.geod.ub;
    
    % Initialize optimizer
    nlp = struct('x', y, 'p', par, 'f', objective, 'g', [ceq;cineq]);
    obj.solver = nlpsol('solver', 'ipopt', nlp, solver_options);
end

function ddelta = f_delta_cont(delta,rho,L_D,w_bar)
    ddelta = -(rho-L_D)*delta + w_bar;
end