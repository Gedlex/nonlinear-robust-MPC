function rho = getOptimalContractionRate(obj,rho,tightenings)
    arguments
        obj (1,1) CCM
        rho (:,:) {mustBeNumeric,mustBeVector}
        tightenings (:,:)
    end
    % Remove nan entries
    idxs = any(isnan(tightenings));
    tightenings(:,idxs) = [];
    rho(idxs) = [];

    % Check terminal constraint
    z_ref = obj.params.x_ref;
    v_ref = obj.params.u_ref(z_ref,zeros(obj.np,1));
    zv_ref = [z_ref; v_ref];

    % Concatenate state and input constraint sets
    A = blkdiag(obj.params.F_x,obj.params.F_u);
    b = [obj.params.b_x; obj.params.b_u];
    
    % Normalize tightenings by terminal constraint
    tightenings = tightenings ./ (b - A*zv_ref);
    
    % Compute sum of normalized tightenings
    norm_sum = sum(tightenings);
    
    % Find minimum
    [~,idx] = min(norm_sum);
    rho = rho(idx);
end