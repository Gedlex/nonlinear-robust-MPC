function obj = compute_delta_over(obj,x_sym,u_sym,theta_sym,d_sym,W_chol_fcn)
    % Determine states used / needed in w_tilde
    w_fcn = obj.w_tilde(x_sym,u_sym,0,theta_sym,d_sym,W_chol_fcn(x_sym));
    vars_used = findDependencies(w_fcn,x_sym,u_sym,theta_sym,d_sym);
    vars_used{1}(obj.W_idxs) = false;
    
    % Compute n_vars
    n_vars = cellfun(@sum, vars_used);
    n_vars = [sum(obj.W_idxs); n_vars(:)];
    
    % Compute number of grid points according to grid pattern
    n_grid = getNgrid(n_vars,obj.grid_pattern,obj.n_tot.delta_over,n_grid=[0,0,~obj.params.param_uncertainty,2,2]);
    
    % Get samples from state and input polyhedron (needed for delta)
    x_param = getSamples(obj.params.F_x,obj.params.b_x,n_grid(1),obj.W_idxs);
    x_eval = getSamples(obj.params.F_x,obj.params.b_x,n_grid(2),vars_used{1});
    u_eval = getSamples(obj.params.F_u,obj.params.b_u,n_grid(3),vars_used{2});
    n_param = size(x_param,2);
    
    % Get all vertex combinations of d and theta (enter delta linearly)
    theta_eval =  getCombinations(obj.params.theta_v.*vars_used{3});
    d_eval = getCombinations(obj.params.disturbance_v.*vars_used{4});
    
    % Create additional variables to prevent communication overhead for parfor
    rho = obj.rho;
    L_D = obj.L_D;
    L_G = obj.L_G;
    W_idxs = obj.W_idxs;
    W_chol_fcn = @(x) W_chol_fcn(x);
    w_tilde = @(x,u,theta,d,W_chol_x) obj.w_tilde(x,u,0,theta,d,W_chol_x);
    
    % Create dataqueue for parfor progbar
    q = parallel.pool.DataQueue;
    afterEach(q, @(~) progbar(increment=100/n_param));
    
    % Initialize progress bar
    progbar(0,prefix='   Delta_over: ')
    
    % Loop over parameter samples
    delta_over = -inf;
    parfor k = 1:n_param
        % Get current state
        x = x_param(:,k);
        
        % Evaluate functions
        W_chol_x = full(W_chol_fcn(x));
        
        % Loop over state samples
        for i = 1:size(x_eval,2)
            % Permute x_eval with x
            x_i = x_eval(:,i);
            x_i(W_idxs) = x(W_idxs);
            
            % Loop over input samples
            for j = 1:size(u_eval,2)
                u_j = u_eval(:,j);
                
                % Loop over theta
                for v = 1:size(theta_eval,2)
                    theta_v = theta_eval(:,v);
                    
                    % Loop over d
                    for l = 1:size(d_eval,2)
                        d_l = d_eval(:,l);
                        
                        % Sample w_tilde
                        delta = w_tilde(x_i,u_j,theta_v,d_l,W_chol_x)/(rho - L_D - L_G*abs(theta_v));
                        
                        % Find maximum delta_over
                        delta_over = max(delta_over, delta);
                    end
                end
            end
        end
        % Update progress bar
        send(q,k);
    end
    obj.delta_over = delta_over;
end