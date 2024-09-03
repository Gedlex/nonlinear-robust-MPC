function y0 = initializeSolver(obj,terminal_constraint,solver_options)
    % Setup nominal MPC
    [solver,bx,bg] = setupMPC(obj,terminal_constraint,solver_options);

    % Define initial solver state for state and input
    z0 = obj.params.x0 + (obj.params.x_ref - obj.params.x0)*linspace(0,1,obj.N+1);
    v0 = obj.params.u_ref(z0,zeros(obj.np,1));
    
    % Initialize solver
    y0 = [repmat(v0,obj.N,1);
          reshape(z0,obj.nx*(obj.N+1),1)];
    
    % Solve optimization
    res = solver('p',obj.params.x0,'x0',y0,'lbx',bx.ylb,'ubx',bx.yub,'lbg',bg.clb,'ubg',bg.cub);
    feasible = solver.stats().success;
    
    if ~feasible
        warning('Solver initialization failed.');
    end
    
    % Extract solution
    y_sol = full(res.x);
    v_sol = y_sol(1:obj.n_v);
    z_sol = y_sol(obj.n_v + (1:obj.n_z));
    
    % Initialize RAMPC solver with nominal trajectory
    idx = obj.n_geod;
    y0 = zeros(obj.n_var,1);
    y0(idx + (1:obj.n_v)) = v_sol;      idx = idx+obj.n_v;
    y0(idx + (1:obj.n_z)) = z_sol;
    if obj.adaptive
        y0(end)= min(obj.delta_bar_x,obj.delta_bar_u);
    end
end