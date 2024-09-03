function obj = compute_c_obs(obj,x_sym,M_chol_fcn,hs_obs_fcn)
    % Determine states used / needed in h_obs
    vars_used = findDependencies(hs_obs_fcn(x_sym),x_sym);
    vars_used{1}(obj.W_idxs) = false;
    
    % Compute n_vars
    n_vars = cellfun(@sum, vars_used);
    n_vars = [sum(obj.W_idxs); n_vars(:)];
    
    % Compute number of grid points according to grid pattern
    n_grid = getNgrid(n_vars,obj.grid_pattern(1:2),obj.n_tot.c_obs);
    
    % Get samples from state polyhedron (needed for h_obs)
    x_param = getSamples(obj.params.F_x,obj.params.b_x,n_grid(1),obj.W_idxs);
    x_eval = getSamples(obj.params.F_x,obj.params.b_x,n_grid(2),vars_used{1});
    n_param = size(x_param,2);

    % Create additional variables to prevent communication overhead for parfor
    W_idxs = obj.W_idxs;
    M_chol_fcn = @(x) M_chol_fcn(x);
    hs_obs_fcn = @(x) hs_obs_fcn(x);

    % Create dataqueue for parfor progbar
    q = parallel.pool.DataQueue;
    afterEach(q, @(~) progbar(increment=100/n_param));
    
    % Initialize progress bar
    progbar(0,prefix='    Obstacles: ')
    
    % Loop over parameter samples
    n_obs = numel(obj.params.h_obs(zeros(obj.nx,1)));
    c_obs = -inf(n_obs,1);
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
            
            % Compute dh_obs/dx
            hs_obs = full(hs_obs_fcn(x_i));
            
            % Find maximum tightening
            c_obs = max(c_obs, vecnorm(hs_obs/M_chol,2,2));
        end
        % Update progress bar
        send(q,k);
    end
    obj.c_obs = c_obs;
end