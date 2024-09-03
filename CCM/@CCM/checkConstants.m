function [d_delta, delta_inf] = checkConstants(obj,z_s0,v_s0,delta_bar_0,M_chol_fcn)
    % In the following, we check if a terminal equality constraint saitsfies the
    % posed conditions with the considered uncertainty. For this we compute the
    % derivative \dot{\delta} at steady-state (for terminal condition).
    
    % Get all vertex combinations of d and theta (enter dynamics linearly)
    theta_eval = getCombinations(obj.params.theta_v);
    d_eval = getCombinations(obj.params.disturbance_v);

    % Create additional variables to prevent communication overhead for parfor
    M_chol_fcn = @(x) M_chol_fcn(x);
    E_fcn = @(x) obj.sys.E_fcn(x);
    G_fcn = @(x,u) obj.sys.G_fcn(x,u);
    
    % Loop over theta
    w_bar_terminal = -inf;
    parfor i = 1:size(theta_eval,2)
        theta_i = theta_eval(:,i);
        
        % Loop over d
        for j = 1:size(d_eval,2)
            d_j = d_eval(:,j);
            
            % Sample w_bar_terminal
            sample = norm(full(M_chol_fcn(z_s0))*(E_fcn(z_s0)*d_j + G_fcn(z_s0,v_s0)*theta_i));
            
            % Find maximum w_bar_terminal
            w_bar_terminal = max(w_bar_terminal,sample);
        end
    end
    
    % Compute d_delta
    theta_max = max(abs(obj.params.theta_v),[],2);
    d_delta = -(obj.rho - obj.L_D - obj.L_G*theta_max)*delta_bar_0 + w_bar_terminal;
    
    if obj.rho_des > 0
        % Compute asymptotic tightening (for bounded tubes)
        delta_inf = w_bar_terminal/(obj.rho - obj.L_D - obj.L_G*theta_max);
        if delta_inf <= 0
            delta_inf = inf;
        end
    else
        % Compute tightening for horizon T (for unbounded tubes)
        T = obj.N*obj.params.dt;
        alpha = (obj.rho - obj.L_D - obj.L_G*theta_max);
        delta_inf = (1 - exp(-alpha*T)) * w_bar_terminal / alpha;
    end
end