function objective = getObjective(obj,z,v,z_ref,v_ref)
    objective = casadi.MX(0);
    for i=1:obj.N
        % Stage cost
        objective = objective + (z(:,i)-z_ref)'*obj.params.Q_cost*(z(:,i)-z_ref) + (v(:,i)-v_ref)'*obj.params.R_cost*(v(:,i)-v_ref);
    end
    
    % Terminal cost
    objective = objective + (z(:,1+obj.N)-z_ref)'*obj.params.P_cost*(z(:,1+obj.N)-z_ref);
end