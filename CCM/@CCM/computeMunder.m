function obj = computeMunder(obj)
    % Compute under-approximation of M that minimizes c_j on position states 
    % for visulaizing an ellipsoidal over-approximation of the tube
    import casadi.*
    
    % Define optimization variables
    yalmip('clear');
    M_under = sdpvar(obj.nx,obj.nx);
    c_j = sdpvar(1);
    
    % Define objective (logdet is nonconvex)
    obj_over = c_j - trace(M_under)*1e-5;
    
    % Define constraints
    con_over = [M_under >= 0]; %#ok
    eps = 1e-10;
    
    % Check if obstacle constraints h_obs are defined
    x_sym = SX.sym('x',obj.nx);
    if isa(obj.params.h_obs,"function_handle")
        % Compute jacobian of h_obs
        hs_obs_fcn = Function('hs_ob',{x_sym},{jacobian(obj.params.h_obs(x_sym),x_sym)});
        
        % Devide total number of sampling points
        n_tot1 = round(obj.n_tot.M_under/2);
        n_tot2 = round(obj.n_tot.M_under/2);
    else
        % Set jacobian to zero
        hs_obs_fcn = Function('hs_ob',{x_sym},{SX.zeros(1,obj.nx)});
        
        % Devide total number of sampling points
        n_tot1 = obj.n_tot.M_under-1;
        n_tot2 = 1;
    end
    
    % Compute n_vars
    n_vars = sum(obj.W_idxs);
    
    % Compute number of grid points according to grid pattern
    n_grid = getNgrid(n_vars,obj.grid_pattern(1),n_tot1);
    
    % Get samples from state polyhedron for parametrization
    x_param = getSamples(obj.params.F_x,obj.params.b_x,n_grid(1),obj.W_idxs);
    
    % Initialize progress bar
    progbar(0,prefix=' Constraint 1: ')
    
    % Loop over parameter samples
    n_param = size(x_param,2);
    for k = 1:n_param
        % Get current state
        x = x_param(:,k);
        
        % Compute inv(W(x))
        M = inv(obj.W_fcn(x));
        
        % M - M_under is positive definite
        con_over = [con_over; M - M_under >= eye(obj.nx)*eps];
        
        % Update progress bar
        progbar(100*k/n_param);
    end
    
    %% Minimizing c_j with M_under    
    % Determine states used / needed in h_obs
    vars_used = findDependencies(hs_obs_fcn(x_sym),x_sym);
    
    % Compute n_vars
    n_vars = cellfun(@sum, vars_used);
    
    % Compute number of grid points according to grid pattern
    n_grid = getNgrid(n_vars,obj.grid_pattern(2),n_tot2);
    
    % Get samples from state polyhedron (needed for h_obs)
    x_eval = getSamples(obj.params.F_x,obj.params.b_x,n_grid(1),vars_used{1});
    
    % Initialize progress bar
    progbar(0,prefix=' Constraint 2: ')
    
    % Loop over state samples
    n_eval = size(x_eval,2);
    for k = 1:n_eval
        % Get current state
        x_i = x_eval(:,k);
        
        % Compute dh_obs/dx
        hs_obs = full(hs_obs_fcn(x_i));
        
        for i = 1:size(hs_obs,1)
            % Construct LMI and add to constraints
            LMI = [M_under,    hs_obs(i,:)';
                   hs_obs(i,:),        c_j];
            con_over = [con_over; LMI >= eye(7)*eps];
        end
        
        % Update progress bar
        progbar(100*k/n_eval);
    end
    
    %% Solve Optimization Problem
    % Solver settings
    ops = sdpsettings('solver','mosek','verbose',0);
    fprintf('Problem formulation finished! Start solving...\n');
    
    % Solve problem
    optimize(con_over,obj_over,ops);
    obj.M_under = value(M_under);
    
    %% Check Obtained Result
    % Compute n_vars
    n_vars = sum(obj.W_idxs);
    
    % Compute number of grid points according to grid pattern
    n_grid = getNgrid(n_vars,obj.grid_pattern(1),obj.n_tot.M_under_check);
    
    % Get samples from state polyhedron for parametrization
    x_param = getSamples(obj.params.F_x,obj.params.b_x,n_grid(1),obj.W_idxs);
    n_param = size(x_param,2);

    % Create additional variables to prevent communication overhead for parfor
    W_fcn = @(x) obj.W_fcn(x);
    M_under = obj.M_under;

    % Create dataqueue for parfor progbar
    q = parallel.pool.DataQueue;
    afterEach(q, @(~) progbar(increment=100/n_param));
    
    % Initialize progress bar
    progbar(0,prefix='Check constraints: ')
    
    % Loop over parameter samples
    M_under_test = inf;
    parfor k = 1:n_param
        % Get current state
        x = x_param(:,k);
        
        % Compute inv(W(x))
        M = inv(W_fcn(x));
        
        % M_under_test should be <= 0
        M_under_test = min(M_under_test, min(eig(M-M_under)));
        
        % Update progress bar
        send(q,k);
    end
    
    if M_under_test > 1e-5
        fprintf(2,['Warning: M_under is inaccurate. Test evaluated to %f.4 (should be < 0). ' ...
                   'Increase number of samples n_tot.M_under which is currently set to %d\n'], ...
                   M_under_test, obj.n_tot.M_under);
    end
end