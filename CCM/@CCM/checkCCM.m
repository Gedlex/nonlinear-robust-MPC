function rho = checkCCM(obj)
    % Check if the contraction rate rho satisfies given CCM conditions.
    
    % Determine variables used / needed in differential dynamics
    dep_A_diff = findDependencies(obj.sys.A_diff,obj.nx,obj.nu,obj.np,obj.nw);
    dep_B_diff = findDependencies(obj.sys.B_diff,obj.nx,obj.nu,obj.np,obj.nw);
    dep_dW_dt  = findDependencies(obj.dW_dt_approx_fcn,obj.nx,obj.nu,obj.np,obj.nw);
    vars_used = cellfun(@(X,Y,Z) X | Y | Z, dep_A_diff, dep_B_diff, dep_dW_dt, 'UniformOutput', false);
    
    % Remove parametrizing states from state samples
    vars_used{1}(obj.W_idxs) = false;
    
    % Compute n_vars
    n_vars = cellfun(@sum, vars_used);
    n_vars = [sum(obj.W_idxs); n_vars(:)];
    
    % Compute number of grid points according to grid pattern
    n_grid = getNgrid(n_vars,obj.grid_pattern,obj.n_tot.rccm_check,n_grid=[0,0,0,2,2]);
    
    % Get samples from state and input polyhedron (needed for differential dynamics)
    x_param = getSamples(obj.params.F_x,obj.params.b_x, n_grid(1), obj.W_idxs);
    x_eval = getSamples(obj.params.F_x,obj.params.b_x, n_grid(2), vars_used{1});
    u_eval = getSamples(obj.params.F_u,obj.params.b_u, n_grid(3), vars_used{2});
        
    % Get all vertex combinations of d and theta (enter dynamics linearly)
    theta_eval = getCombinations(obj.params.theta_v.*vars_used{3});
    d_eval = getCombinations(obj.params.disturbance_v.*vars_used{4});
    n_param = size(x_param,2);

    % Create additional variables to prevent communication overhead for parfor
    W_idxs = obj.W_idxs;
    W_fcn = @(x) obj.W_fcn(x);
    Y_fcn = @(x) obj.Y_fcn(x);
    A_diff = @(x,u,theta,d) obj.sys.A_diff(x,u,theta,d);
    B_diff = @(x,u,theta,d) obj.sys.B_diff(x,u,theta,d);
    dW_dt_approx_fcn = @(x,u,theta,d) obj.dW_dt_approx_fcn(x,u,theta,d);

    % Create dataqueue for parfor progbar
    q = parallel.pool.DataQueue;
    afterEach(q, @(~) progbar(increment=100/n_param));
    
    % Initialize contraction rate and eigenvalues
    rho = inf;
    
    % Initialize progress bar
    progbar(0,prefix='Check CCM: ')
    
    % Loop over parameter samples
    parfor k = 1:n_param
        % Get current state
        x = x_param(:,k);
        
        % Evaluate functions
        W_k = W_fcn(x);
        Y_k = Y_fcn(x);
        
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
                        
                        % Our condition on transformed CCM
                        %-W_dot + <A_cl W> + 2 rho W >= 0 <--> -W_dot + <A W + B_hat Y> + 2 rho W >= 0
                        tmp = A_diff(x_i,u_j,theta_v,d_l)*W_k + B_diff(x_i,u_j,theta_v,d_l)*Y_k;
                        F1 = -dW_dt_approx_fcn(x_i,u_j,theta_v,d_l) + (tmp + tmp');
                        
                        % Set rho based on generalized eigenvalue, i.e., biggest rho such that F_1+2*\rho*W<=0
                        rho = min(rho, -max(eig(F1,W_k))/2); % maximal ev. should be <= 0
                    end
                end
            end
        end
        % Update progress bar
        send(q,k);
    end
   fprintf('Contraction rate: %.4f\n', rho);
end