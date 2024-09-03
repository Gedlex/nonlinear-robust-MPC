function obj = compute_L_G(obj,x_sym,u_sym,theta_sym,M_chol_fcn,Gs_fcn)
    % Determine states used / needed in Gs
    vars_used = findDependencies(Gs_fcn(x_sym,u_sym,theta_sym),x_sym,u_sym,theta_sym);
    vars_used{1}(obj.W_idxs) = false;
    
    % Compute n_vars
    n_vars = cellfun(@sum, vars_used);
    n_vars = [sum(obj.W_idxs); n_vars(:)];
    
    % Compute number of grid points according to grid pattern
    n_grid = getNgrid(n_vars,obj.grid_pattern(1:4),obj.n_tot.L_G,n_grid=[0,0,0,1]);
    
    % Get samples from state and input polyhedron (needed for L_G)
    x_param = getSamples(obj.params.F_x,obj.params.b_x,n_grid(1),obj.W_idxs);
    x_eval = getSamples(obj.params.F_x,obj.params.b_x,n_grid(2),vars_used{1});
    u_eval = getSamples(obj.params.F_u,obj.params.b_u,n_grid(3),vars_used{2});
    n_param = size(x_param,2);

    % Create additional variables to prevent communication overhead for parfor
    W_idxs = obj.W_idxs;
    M_chol_fcn = @(x) M_chol_fcn(x);
    Gs_fcn = @(x,u,theta_idx) Gs_fcn(x,u,theta_idx);
    np = obj.np;
    
    % Create dataqueue for parfor progbar
    q = parallel.pool.DataQueue;
    afterEach(q, @(~) progbar(increment=100/n_param));
    
    % Initialize progress bar
    progbar(0,prefix='          L_G: ')
    
    % Loop over parameter samples
    L_G = -inf(obj.np,n_param);
    parfor k = 1:n_param
        % Get current state
        x = x_param(:,k);
        
        % Compute Cholesky decomposition of M(x)
        M_chol = full(M_chol_fcn(x));
        
        % Loop over state samples
        for i = 1:size(x_eval,2)
            % Permute x_eval with x
            x_i = x_eval(:,i);
            x_i(W_idxs) = x(W_idxs);
            
            % Loop over input samples
            for j = 1:size(u_eval,2)
                u_j = u_eval(:,j);
                
                % Loop over all theta
                for theta_idx = 1:np
                    theta_i = false(np,1);
                    theta_i(theta_idx) = true;
                    
                    % Compute Gs
                    Gs = full(Gs_fcn(x_i,u_j,theta_i));
                    
                    % Find maximum continuity constant
                    L_G(theta_idx,k) = max(L_G(theta_idx,k), norm(Gs/M_chol,2));
                end
            end
        end
        % Update progress bar
        send(q,k);
    end
    obj.L_G = max(L_G,[],2);
end