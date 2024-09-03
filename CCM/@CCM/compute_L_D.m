function obj = compute_L_D(obj,x_sym,d_sym,M_chol_fcn,Es_fcn)
    % Determine states used / needed in Es
    vars_used = findDependencies(Es_fcn(x_sym,d_sym),x_sym,d_sym);
    vars_used{1}(obj.W_idxs) = false;
    
    % Compute n_vars
    n_vars = cellfun(@sum, vars_used);
    n_vars = [sum(obj.W_idxs); n_vars(:)];
    
    % Compute number of grid points according to grid pattern
    n_grid = getNgrid(n_vars,obj.grid_pattern([1,2,5]),obj.n_tot.L_D,n_grid=[0,0,2]);
    
    % Get samples from state polyhedron (needed for L_D)
    x_param = getSamples(obj.params.F_x,obj.params.b_x,n_grid(1),obj.W_idxs);
    x_eval = getSamples(obj.params.F_x,obj.params.b_x,n_grid(2),vars_used{1});
    n_param = size(x_param,2);
    
    % Get all vertex combinations of d (enter L_D linearly)
    d_eval = getCombinations(obj.params.disturbance_v.*vars_used{2});

    % Create additional variables to prevent communication overhead for parfor
    W_idxs = obj.W_idxs;
    M_chol_fcn = @(x) M_chol_fcn(x);
    Es_fcn = @(x,d) Es_fcn(x,d);

    % Create dataqueue for parfor progbar
    q = parallel.pool.DataQueue;
    afterEach(q, @(~) progbar(increment=100/n_param));
    
    % Initialize progress bar
    progbar(0,prefix='          L_D: ')
    
    % Loop over parameter samples
    L_D = -inf;
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
    
            % Loop over d
            for j = 1:size(d_eval,2)
                d_j = d_eval(:,j);
    
                % Compute Es
                Es = full(Es_fcn(x_i,d_j));
                
                % Find maximum continuity constant
                L_D = max(L_D, norm(Es/M_chol,2));
            end
        end
        % Update progress bar
        send(q,k);
    end
    obj.L_D = L_D;
end