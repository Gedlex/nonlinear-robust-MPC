function delta_tight = deltaProp(obj,z,v,delta,theta_bar,theta_v)
    % Propagate tube scaling
    % Tube scaling delta in MPC can also be an over-approximation, e.g., if the
    % constraints are not active the optimal tube scaling is not unique 
    % -> for visualization we re-compute the smallest scaling that satisfes the conditions with equality.
    
    % Define joint dynamics of s = [z, delta];
    f = @(s,v) f_joint(obj,s,v,theta_bar,theta_v);
    
    % Allocate joint state s
    s = zeros(obj.nx+1,obj.N+1);
    s(:,1) = [z(:,1);delta(1)];
    
    % Propagate delta
    for i = 2:obj.N+1
        s(:,i) = obj.sys.ddyn(s(:,i-1),v(:,i-1),dynamics=f);
        s(1:obj.nx,i) = z(:,i);
    end
    
    % Extract delta tight
    delta_tight = s(end,:)';
end

function dot_joint = f_joint(obj,s,v,theta_bar,theta_v)
    % Extract joint states
    z = s(1:obj.nx);
    delta = s(obj.nx+1);
    
    % Get all vertex combinations of d and theta (enter delta dynamics linearly)
    theta_eval =  getCombinations(theta_v);
    d_eval = getCombinations(obj.params.disturbance_v);
    
    % Loop over theta
    max_mismatch = -inf;
    for i = 1:size(theta_eval,2)
        theta_i = theta_eval(:,i);
        
        % Loop over d
        for j = 1:size(d_eval,2)
            d_j = d_eval(:,j);
            
            % Compute w_tilde
            weighted_norm = obj.w_tilde(z,v,theta_bar,theta_i,d_j,chol(obj.W_fcn(z)));

            % Find maximum mismatch
            max_mismatch = max(max_mismatch, obj.L_G*abs(theta_i - theta_bar)*delta + weighted_norm);
        end
    end
    
    % Dynamics for z and delta
    dot_joint(1:obj.nx,1) = obj.sys.fw(z,v,theta_bar,0);
    dot_joint(obj.nx+1,1) = -(obj.rho - obj.L_D)*delta + max_mismatch;
end