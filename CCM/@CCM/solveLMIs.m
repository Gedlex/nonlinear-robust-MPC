function [obj, success] = solveLMIs(obj,rho,v_W_fun,dv_W_dx_fun,W_coef,Y_coef)
    
    % Determine variables used / needed in LMIs
    v_fw = @(x,u,theta,d) double(obj.state2param)*obj.sys.fw(x,u,theta,d);
    dep_v_fw = findDependencies(v_fw,obj.nx,obj.nu,obj.np,obj.nw);
    dep_A_diff = findDependencies(obj.sys.A_diff,obj.nx,obj.nu,obj.np,obj.nw);
    dep_B_diff = findDependencies(obj.sys.B_diff,obj.nx,obj.nu,obj.np,obj.nw);
    dep_BwXw = findDependencies(obj.sys.BwXw,obj.nx,obj.nu,obj.np,obj.nw);
    
    vars_used = cellfun(@(V,X,Y,Z) V | X | Y | Z, dep_v_fw, dep_A_diff, dep_B_diff, dep_BwXw, 'UniformOutput', false);
    
    % Remove parametrizing states from state samples
    vars_used{1}(obj.W_idxs) = false;
    
    % Compute n_vars = [n_params; n_x_eval; n_u_eval; n_theta_eval; n_d_eval]
    n_vars = cellfun(@sum, vars_used);
    n_vars = [sum(obj.W_idxs); n_vars(:)];
    
    % Compute number of grid points according to grid pattern
    n_grid = getNgrid(n_vars,obj.grid_pattern,round(obj.n_tot.rccm/3),n_grid=[0,0,0,2,2]);
    
    % Get samples from state and input polyhedron (needed for differential dynamics)
    x_param = getSamples(obj.params.F_x,obj.params.b_x, n_grid(1), obj.W_idxs);
    x_eval = getSamples(obj.params.F_x,obj.params.b_x, n_grid(2), vars_used{1});
    u_eval = getSamples(obj.params.F_u,obj.params.b_u, n_grid(3), vars_used{2});
    
    % Get all vertex combinations of d and theta (enter dynamics linearly)
    theta_eval = getCombinations(obj.params.theta_v.*vars_used{3});
    d_eval = getCombinations(obj.params.disturbance_v.*vars_used{4});
    n_param = size(x_param,2);

    % Define hyperparameters
    lambda = 2*rho;         % From Zhao's paper
    scaling = 10^(-2.5);    % Cost scaling
    yScaling = 2;           % Scaling for Y vars
    wScaling = 1;           % Scaling for W vars
    wbar_lb = 1E-6;         % Lower bound for wbar
    W_lb = 1e-2;            % Lower bound for W
    
    % Create optimization variables
    wbar = sdpvar;
    W_lower = sdpvar;

    % Define objective
    vars_wd = [W_coef(:)*wScaling; Y_coef(:)*yScaling];
    objective = norm(vars_wd,1)*scaling + wbar;

    % Define constraints
    eps = 1E-6;
    constraints = [wbar>=wbar_lb;
                   W_lower >= W_lb];

    % Create additional variables to prevent communication overhead for parfor
    W_idxs = obj.W_idxs;
    A_diff = @(x,u,theta,d) obj.sys.A_diff(x,u,theta,d);
    B_diff = @(x,u,theta,d) obj.sys.B_diff(x,u,theta,d);
    BwXw   = @(x,u,theta,d) obj.sys.BwXw(x,u,theta,d);
    nx = obj.nx;
    nu = obj.nu;
    
    % Initialize progress bar
    progbar(0,prefix='Construct LMIs: ')
    
    % Loop over parameter samples
    for k = 1:n_param
        % Get current state
        x = x_param(:,k);

        % Compute monomials and their derivative (w.r.t. x)
        v_W = v_W_fun(x);
        dv_W_dx = dv_W_dx_fun(x);
        n_monos_W = length(v_W);

        % Parametrize Y
        Y = zeros(nu,nx);
        for i = 1:n_monos_W
            Y = Y + Y_coef(:,:,i)*v_W(i);
        end
        
        % Parametrize W
        W = zeros(nx);
        for i = 1:n_monos_W
            W = W + W_coef(:,:,i)*v_W(i);
        end
        
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
                        
                        % Compute time derivative dv_W_dt
                        dv_W_dt = dv_W_dx*v_fw(x_i,u_j,theta_v,d_l);
                        
                        % Parametrize dW_dt
                        dW_dt = zeros(nx);
                        for m = 1:n_monos_W
                            dW_dt = dW_dt + W_coef(:,:,m)*dv_W_dt(m); 
                        end
                        
                        % LMI to ensure contraction rate lambda
                        ddyn = A_diff(x_i,u_j,theta_v,d_l)*W + B_diff(x_i,u_j,theta_v,d_l)*Y;
                        LMI1 = dW_dt - (ddyn+ddyn') - lambda*W;
                        
                        % LMI to ensure |BwXw|_M <= wbar
                        BwXw_l = BwXw(x_i,u_j,theta_v,d_l);
                        LMI2 = [W,       BwXw_l; ...
                                BwXw_l',  wbar];
                        
                        % LMI to ensure that W_lower <= W
                        LMI3 =  W - W_lower*eye(nx);
                        
                        % Add constraints
                        constraints = [constraints;
                                       LMI1 >= eye(nx)*eps;
                                       LMI2 >= eye(nx+1)*eps;
                                       LMI3 >= eye(nx)*eps];
                    end
                end
            end
        end
        % Update progress bar
        progbar(100*k/n_param);
    end
    
    %% Solve Optimization Problem
    % Solver settings
    ops = sdpsettings('solver','mosek','verbose',0);
    fprintf('LMI formulation finished! Start solving...\n');
    
    % Solve problem
    sol = optimize(constraints,objective,ops);
    
    % Check solution
    if (sol.problem == 0 || sol.problem == 4) 
        wbar_opt = value(wbar);
        success = true;
    else
        wbar_opt = inf;
        success = false;
    end
    
    % Display cost and wbar_opt
    vars_wd = [value(W_coef(:))*wScaling; value(Y_coef(:))*yScaling];
    fprintf('Coef Cost:  %.3f\nw_bar Cost: %.5f\n', norm(vars_wd,1)*scaling, wbar_opt);
    fprintf('RCCM rho = %.4f, w_bar = %.4f\n', rho, wbar_opt);
end
