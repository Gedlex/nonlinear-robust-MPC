function y0 = shiftSolution(obj,y_sol,delta_bar_sol,theta_bar_sol)
    % Remain on current trajectory by shifting solution.
    
    % Initialize shifted y0
    y0 = y_sol;
    
    % Get steady-state input
    u_eq = obj.params.u_ref(obj.params.x_ref,theta_bar_sol);
    
    % Shift inputs
    idx = obj.n_geod;
    y0(idx + (1:obj.n_v-obj.nu)) = y_sol(idx + (obj.nu+1:obj.n_v)); idx = idx+obj.n_v-obj.nu;
    y0(idx + (1:obj.nu)) = u_eq;                                    idx = idx+obj.nu;
    
    % Shift states
    y0(idx + (1:obj.n_z-obj.nx)) = y_sol(idx + (obj.nx+1:obj.n_z)); idx = idx+obj.n_z-obj.nx;
    y0(idx + (1:obj.nx)) = y_sol(idx + (1:obj.nx));                 idx = idx+obj.nx;
    
    % Shift tube scaling
    y0(idx + (1:obj.n_delta-1)) = y_sol(idx + (1+1:obj.n_delta));   idx = idx+obj.n_delta-1;
    y0(idx+1) = delta_bar_sol;                                      idx = idx+1;
    
    % Shift w_bar (4 because of rk4)
    y0(idx + (1:obj.n_w_bar-4)) = y_sol(idx + (4+1:obj.n_w_bar));   idx = idx+obj.n_w_bar-4;
    y0(idx + (1:4)) = y_sol(idx + (1:4));
end
