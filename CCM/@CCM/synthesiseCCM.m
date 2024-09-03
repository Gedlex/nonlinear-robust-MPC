function [obj, success] = synthesiseCCM(obj,rho)
    %% Parametize Y, W and dW_dt
    % Create symbolic variables and extract parametrization
    x = sym('x', [obj.nx,1]);
    u = sym('u', [obj.nu,1]);
    d = sym('d', [obj.nw,1]);
    theta = sym('theta',[obj.np,1]);
    W_states = obj.state2param*x;
    
    % Get monomials v_W and their derivative dv_W_dx
    v_W = sym(monolist(W_states, obj.monomial_degree));
    dv_W_dx = jacobian(v_W, W_states);
    n_monos_W = length(v_W);
    
    % Convert to matlab function
    v_W_fun = matlabFunction(v_W,'Vars',{x});
    dv_W_dx_fun = matlabFunction(dv_W_dx,'Vars',{x});
    
    % Parametrize Y and W
    yalmip('clear');
    Y_coef = sdpvar(obj.nu,obj.nx,n_monos_W);
    W_coef = sdpvar(obj.nx,obj.nx,n_monos_W);
    
    %% Construct CCM
    tic;
    if obj.use_sos
        [obj, success] = obj.solveSOS(rho,v_W_fun,dv_W_dx_fun,W_coef,Y_coef);
    else
        [obj, success] = obj.solveLMIs(rho,v_W_fun,dv_W_dx_fun,W_coef,Y_coef);
    end
    if ~success; return; end
    toc;
    
    %% Extract Solution
    % Get coefficients
    Y_coef = value(Y_coef);
    obj.Y_coef = clean(Y_coef,1e-10);
    
    % Extract Y_fcn
    Y = zeros(obj.nu,obj.nx);
    for i=1:n_monos_W
        Y = Y + obj.Y_coef(:,:,i)*v_W(i);
    end
    
    % Convert to matlab function
    obj.Y_fcn = matlabFunction(Y,'Vars',{x});
    
    % Get coefficients
    W_coef = value(W_coef);
    obj.W_coef = clean(W_coef,1e-10);
    
    % Extract W_fcn
    W = zeros(obj.nx);
    for i=1:n_monos_W
        W = W + obj.W_coef(:,:,i)*v_W(i);
    end
    
    % Convert to matlab function
    obj.W_fcn = matlabFunction(W,'Vars',{x});  
    
    % Extract dW_dx
    m = sum(obj.W_idxs);
    idxs = find(obj.W_idxs == 1);
    dW_dx = repmat({zeros(obj.nx)},obj.nx,1);
    for i = 1:m
        idx = idxs(i);
        for j=1:n_monos_W
            dW_dx{idx} = dW_dx{idx} + obj.W_coef(:,:,j)*dv_W_dx(j,i); 
        end
    end
    
    % Extract dynamics for parametrization
    w_dyn = obj.state2param*obj.sys.fw(x,u,theta,d);
    w_dyn_approx = obj.state2param*obj.sys.fw_approx(x,u,theta,d);
    
    % Compute dW_dt (exact and approximate)
    dW_dt = zeros(obj.nx,obj.nx);
    dW_dt_approx = zeros(obj.nx,obj.nx);
    for i = 1:m
        idx = idxs(i);
        dW_dt = dW_dt + dW_dx{idx}*w_dyn(i);
        dW_dt_approx = dW_dt_approx + dW_dx{idx}*w_dyn_approx(i);
    end
    
    % Convert to matlab function
    obj.dW_dx = matlabFunction(dW_dx{:},'Vars',{x});
    obj.dW_dt_fcn = matlabFunction(dW_dt,'Vars',{x,u,theta,d});
    obj.dW_dt_approx_fcn = matlabFunction(dW_dt_approx,'Vars',{x,u,theta,d});
    
    %% Check CCM Condition
    startParpool();
    rho_eval = obj.checkCCM();
    
    % Return success or fail
    if rho_eval <= 0 && rho > 0
        success = false;
    else
        obj.rho = rho_eval;
        success = true;
        if rho <= 0
            fprintf(2,'Warning: The rho obtained is nonpositive, hence the tube is unbounded. Check your solution carefully.')
        end
    end
end