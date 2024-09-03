function [obj, success] = solveSOS(obj,rho,v_W_fun,dv_W_dx_fun,W_coef,Y_coef)
    %  Compute W and Y
    %
    %  Based on the code of Pan Zhao, UIUC, Advanced Controls Research Lab,
    %  panzhao2@illinois.edu
    %  https://github.com/boranzhao/robust_ccm_tube
    %  for the paper:
    %  P. Zhao, et al. Tube-certified trajectory tracking for nonliner systems
    %  with robust control contraction metrics. Submitted to IEEE Robotics and
    %  Automation Letters, 2021
    %  Last update: Sep 9, 2021
    %
    %  Compute W and Y satisfying condition (15) in paper using SOS programming  
    %  with S procedure.
    %  
    %  Changes: Account for parametric uncertainty (theta) as well
    %  - B changed to B_hat (dependent on theta)
    %  - B_w depends on u
    %  - Condition (15) depends on u and theta now -> modified S procedure
    
    % x is the variable that appear in W and Y
    
    warning('solveSOS is legacy code which is only compatible with the Planar Quadrotor.')
    
    % Define hyperparameters
    tolerance = 1e-6;         % Tolerance or wbar
    lambda = 2*rho;           % From Zhao's paper
    W_lower_bound = 1e-2;     % Enforce W>=1e-2
    lagrange_deg_W = 4;       % Bound of W
    lagrange_deg_ccm = 4;     % Ccm condition; had to reduce since there are 3 new variables: theta, u(1) and u(2)
    yScaling = 2;             % Scaling for Y vars
    wScaling = 1;             % Scaling for W vars
    scaling = 10^(-2.5);      % Cost scaling
    
    % Create optimization variables
    x = sdpvar(obj.nx,1);     % State
    d = sdpvar(obj.nw,1);     % Disturbance
    theta = sdpvar(obj.np,1); % Parametric uncertainty
    u = sdpvar(obj.nu,1);     % Input
    y = sdpvar(obj.nx+1,1);   % Variables to convert matrix-valuded SOS to scalar-valued SOS
    wbar = sdpvar;            % Disturbance upper bound
    W_lower = sdpvar;         % W lower bound
    
    % Extract variables
    y1 = y(1:obj.nx);
    yW = y(1:obj.nx);
    W_states = obj.state2param*x;
    
    % Define parameters and constraints
    constraints = [];
    params = [W_coef(:);Y_coef(:);W_lower];

    % Get monomials
    v_W = v_W_fun(x);
    dv_W_dx = dv_W_dx_fun(x);
    n_monos_W = length(v_W);
    
    % Parametrize Y
    Y = zeros(obj.nu,obj.nx);
    for i=1:n_monos_W
        Y = Y + Y_coef(:,:,i)*v_W(i);
    end
    
    % Parametrize W
    W = zeros(obj.nx);
    for i=1:n_monos_W
        W = W + W_coef(:,:,i)*v_W(i);
    end
    
    % Compute time derivative dv_W_dt
    dv_W_dt = dv_W_dx*double(obj.state2param)*obj.sys.fw_approx(x,u,theta,d);  
    
    % Parametrize dW_dt
    dW_dt = zeros(obj.nx);
    for i=1:n_monos_W
        dW_dt = dW_dt + W_coef(:,:,i)*dv_W_dt(i); 
    end
    
    % Get 1D ellipsoidal box constraints
    box_lim_fun = poly2box(obj.params.F_x(5:end,:),obj.params.b_x(5:end), ...% Exclude position constraints
                           obj.params.F_u,obj.params.b_u, ...
                           obj.params.F_d,obj.params.b_d, ...
                           obj.params.theta_v);
    box_lim = box_lim_fun(x,u,theta,d);
    other_lim_states = [x(6);x(5)];
    
    % Get differential dynamics
    A_diff = obj.sys.A_diff(x,u,theta,d);
    B_diff = obj.sys.B_diff(x,u,theta,d);
    BwXw = obj.sys.BwXw(x,u,theta,d);
    
    % Define LMI to ensure contraction rate lambda
    tmp = A_diff*W + B_diff*Y; % B_hat instead of B
    M_pos1 = dW_dt-(tmp+tmp')-lambda*W;
    yMy1 = y1'*M_pos1*y1;
    
    % Define LMI to ensure |BwXw|_M <= wbar (rewritten using Schur complement)
    M_pos2 = [W,     BwXw; ...
              BwXw', wbar];
    yMy2 = y'*M_pos2*y;
    
    % Define LMI to ensure that W_lower <= W
    y_W_Wlower_y =  yW'*(W-W_lower*eye(obj.nx))*yW;
    
    % Initialize progress bar
    progbar(0,prefix='Construct SOS: ')
    
    % Apply s-procedure using compact sets
    n_box = length(box_lim);
    for i=1:n_box
        % W depends on x_3, x_4 -> 
        % Corresponding entries in box_lim: 1,2
        if i==1 || i==2  
            % Enforce a lower bound for W            
            [~,c_lower,v_lower] = polynomial([W_states;yW],lagrange_deg_W);
            
            % Only take the terms quadratic in y
            index = [];
            
            for k=1:length(v_lower)
                if sum(degree(v_lower(k),yW)) == 2
                    index = [index k];
                end
            end
            c_lower = c_lower(index); v_lower = v_lower(index);
            L_lower = c_lower'*v_lower;            
            y_W_Wlower_y =  y_W_Wlower_y-L_lower*box_lim(i);
    
        end
        
        % Contraction condition
        % Depends on x_3, x_4, x_5, x_6 and w_1, w_2 ->
        % Corresponding entries in box_lim: 1,2,4,3,7,8
        if i == 1 || i ==2 || i ==3 || i == 4 || i ==7 || i == 8
            variables = [W_states;other_lim_states;theta;d;y1];
            [~,c_L1,v_L1] = polynomial(variables,lagrange_deg_ccm);
    
            % Only take the terms quadratic in y
            index = [];
            for k=1:length(v_L1)
                if sum(degree(v_L1(k),y1)) == 2
                    index = [index k];
                end
            end
            c_L1 = c_L1(index); v_L1 = v_L1(index);
            
            L1 = c_L1'*v_L1;
            
            yMy1 = yMy1 - L1*box_lim(i);        
            params = [params;vec(c_L1)];
            constraints = [constraints sos(L1)];
        end
        
        % Second LMI depends on x_3, x_4, u_1, u_2 and w_1, w_2 -> 
        % Corresponding entries in box_lim: 1,2,5,6,7,8
        if i == 1 || i ==2 || i ==5 || i == 6 || i ==7 || i == 8
            [~,c_L2,v_L2] = polynomial([W_states;theta;d;u;y],lagrange_deg_W);
    
            % Only take the terms quadratic in y
            index = [];
            for k=1:length(v_L2)
                if sum(degree(v_L2(k),y)) == 2
                    index = [index k];
                end
            end
            c_L2 = c_L2(index); v_L2 = v_L2(index);
            
            L2 = c_L2'*v_L2;
            
            yMy2 = yMy2 - L2*box_lim(i);
            params = [params;vec(c_lower);vec(c_L2)];
            constraints = [constraints sos(L_lower) sos(L2)];
        end
        % Update progress bar
        progbar(100*i/n_box);
    end
    
    % SOS constraints
    constraints = [constraints  sos(yMy1) sos(yMy2) sos(y_W_Wlower_y) wbar>=tolerance W_lower >= W_lower_bound];

    % SOS solver settings
    ops = sdpsettings('solver','mosek','verbose',0);
    fprintf('Problem formulation finished! Start solving...\n');
    
    vars_wd = [W_coef(:)*wScaling; Y_coef(:)*yScaling];
    objective = norm(vars_wd,1) * scaling + wbar;

    % Solve SOS minimizing maximial disturbance bound wbar and coefficients in WY
    [sol,~,~,~] = solvesos(constraints,objective,ops,params);
    disp(sol.info)

    % Check solution
    if (sol.problem == 0 || sol.problem == 4) 
        wbar_opt = value(wbar);
        success = true;
    else
        wbar_opt = inf;
        success = false;
    end
    fprintf('RCCM rho = %.4f, w_bar = %.4f\n', rho, wbar_opt);
end