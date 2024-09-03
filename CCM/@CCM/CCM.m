% -------------------------------------------------------------------------
% File: CCM.m
% Author: Alexander Erdin (aerdin@ethz.ch)
% Date: 19th April 2024
% License: MIT
% Reference:
% Sasfi A.,Köhler J., Zeilinger MN. "Robust adaptive MPC using control
% contraction metrics".
% Link: https://arxiv.org/abs/2209.11713
% -------------------------------------------------------------------------
classdef CCM
    properties (GetAccess=public,SetAccess=private)
        sys (1,1) {mustBeA(sys,["NonlinearSystem","double"])} = NaN        % System dynamics
        params (1,1) {mustBeA(params,["Parameters","double"])} = NaN       % Parameters
        adaptive (1,1) logical                                             % Boolean for set-membership estimaton
        
        N (1,1) {mustBeInteger,mustBeNonnegative}                          % Horizon (prediction steps)
        nx (1,1) {mustBeInteger,mustBeNonnegative}                         % Number of states
        nu (1,1) {mustBeInteger,mustBeNonnegative}                         % Number of inputs
        np (1,1) {mustBeInteger,mustBeNonnegative}                         % Number of parametric uncertainty
        nw (1,1) {mustBeInteger,mustBeNonnegative}                         % Number of disturbance states
        
        W_idxs (:,:) logical {mustBeVector(W_idxs,"allow-all-empties")}    % States parametrizing CCM
        state2param (:,:) logical                                          % Matrix to transform x -> x_param
        monomial_degree (1,1) {mustBeInteger}                              % Monomials parametrizing W and Y
        rho (1,1) {mustBeNumeric} = inf                                    % Contraction rate
        rho_des (1,1) {mustBeNumeric} = inf                                % Desired contraction rate
        use_sos (1,1) logical                                              % Boolean indicating how rccm is constructed

        Y_fcn (1,1) {mustBeA(Y_fcn,["function_handle","double"])} = NaN    % Y
        Y_coef (:,:,:) {mustBeNumeric}                                     % Monomial coefficients of Y
        W_fcn (1,1) {mustBeA(W_fcn,["function_handle","double"])} = NaN    % W
        W_coef (:,:,:) {mustBeNumeric}                                     % Monomial coefficients of W
        dW_dx (1,1) {mustBeA(dW_dx,["function_handle","double"])} = NaN    % State derivative of W
        dW_dt_fcn (1,1) {mustBeA(dW_dt_fcn,["function_handle","double"])} = NaN % Time derivative of W
        dW_dt_approx_fcn (1,1) {mustBeA(dW_dt_approx_fcn,["function_handle","double"])} = NaN % Approximate time derivative of W

        c_x (:,:) {mustBeNumeric,mustBeVector(c_x,"allow-all-empties")}    % Tightening on state constraints
        c_u (:,:) {mustBeNumeric,mustBeVector(c_u,"allow-all-empties")}    % Tightening on input constraints
        c_obs (:,:) {mustBeNumeric,mustBeVector(c_obs,"allow-all-empties")}% Tightening on obstacle constraints
        L_D (1,1) {mustBeNumeric}                                          % Continuity constant for disturbance
        L_G (:,:) {mustBeNumeric,mustBeVector(L_G,"allow-all-empties")}    % Continuity constant for parametric uncertainty
        M_under (:,:) {mustBeNumeric}                                      % Underapproximation for M
        delta_bar_x (:,:) {mustBeNumeric}                                  % Steady-state tightening on state constraints
        delta_bar_u (:,:) {mustBeNumeric}                                  % Steady-state tightening on input constraints
        delta_over (:,:) {mustBeNumeric}                                   % Delta_over for rigid tube mpc
        delta_inf (:,:) {mustBeNumeric}                                    % Asymptotic tightening
        terminal_invariance (1,1) logical                                  % Boolean indicating if terminal invariance was achieved
        terminal_constraint (1,1) logical                                  % Boolean indicating if terminal constraints are used
        exact_initialization (1,1) logical                                 % Boolean indicating if delta is initialized at zero
        
        solver (1,1) {mustBeA(solver,["casadi.Function","double"])} = NaN  % NLP optimizer
    end
    
    properties (Access=public)
        y0 (:,:) {mustBeNumeric,mustBeVector(y0,"allow-all-empties")}      % Initial value of solver
    end
    
    properties (GetAccess=public,SetAccess=private,Hidden)
        n_tot (1,1) struct                                                 % Number of grid points in total
        grid_pattern (5,1) {mustBeNonnegative}                             % Number of grid points relative to each other
        n_eq (1,1) {mustBeInteger,mustBeNonnegative}                       % Number of equalities
        n_ineq (1,1) {mustBeInteger,mustBeNonnegative}                     % Number of inequalities
        n_var (1,1) {mustBeInteger,mustBeNonnegative}                      % Number of optimization variables
        n_v (1,1) {mustBeInteger,mustBeNonnegative}                        % Number of nominal input variables
        n_z (1,1) {mustBeInteger,mustBeNonnegative}                        % Number of nominal state variables
        n_delta (1,1) {mustBeInteger,mustBeNonnegative}                    % Number of contraction variables
        n_w_bar (1,1) {mustBeInteger,mustBeNonnegative}                    % Number of distruabnce variables
        n_geod (1,1) {mustBeInteger,mustBeNonnegative}                     % Number of geodesic variables
        xineq_idxs (:,1) {mustBeInteger,mustBeNonnegative}                 % Indices of state inequality constraints
        uineq_idxs (:,1) {mustBeInteger,mustBeNonnegative}                 % Indices of input inequality constraints
        geod (1,1) struct                                                  % Geodesic struct
        clb (:,:) {mustBeNumeric,mustBeVector(clb,"allow-all-empties")}    % Lower bounds on constraints
        cub (:,:) {mustBeNumeric,mustBeVector(cub,"allow-all-empties")}    % Upper bounds on constraints
        ylb (:,:) {mustBeNumeric,mustBeVector(ylb,"allow-all-empties")}    % Lower bounds on decision variables
        yub (:,:) {mustBeNumeric,mustBeVector(yub,"allow-all-empties")}    % Upper bounds on decision variables
    end
    
    methods
        function obj = CCM(sys,params,options)
            arguments
                sys (1,1) NonlinearSystem
                params (1,1) Parameters
                options.W_idxs (:,:) {mustBeInteger,mustBePositive,mustBeVector} = [3,4]
                options.rho (:,:) {mustBeVector} = 0.7
                options.monomial_degree (1,1) {mustBeInteger,mustBeNonnegative} = 4
                options.use_sos (1,1) logical = false
                options.adaptive (1,1) logical = false
                options.terminal_constraint (1,1) logical = true
                options.exact_initialization (1,1) logical = true
                options.recompute_rccm (1,1) logical = false
                options.recompute_constants (1,1) logical = false;
                options.grid_pattern (5,1) {mustBePositive} = [2,1,1,1,1]
                options.solver_options (1,1) struct = struct('ipopt', struct('max_iter',10000, 'print_level',0), 'print_time',0)
                options.n_tot (1,1) struct = struct('rccm',1E4,'rccm_check',1E5,'c_x',1E5,'c_u',1E5,'c_obs',1E5,'L_D',1E5,'L_G',1E5,'w_bar',1E5,'delta_over',1E5,'M_under',5E3,'M_under_check',1E4)
                options.recompute (1,1) struct = struct('rccm',false,'c_x',false,'c_u',false,'c_obs',false,'L_D',false,'L_G',false,'delta_over',false,'M_under',false)
            end
            % Check for casadi
            if ~exist('casadi.MX', 'class')
                id = 'CCM:UndefinedClass';
                error(id,'CASADI not found. Add CASADI to the MATLAB search path.')
            end
            import casadi.*
            
            % Check for yalmip
            if ~exist('sdpvar','file')
                id = 'CCM:UndefinedClass';
                error(id,'Yalmip not found. Add Yalmip to the MATLAB search path.')
            end
            
            % Check for cprnd
            if ~exist('cprnd', 'file')
                id = 'CCM:UndefinedFunction';
                error(id,'CPRND not found. Add CPRND to the MATLAB search path.')
            end
            
            % Check that sampling intervals match
            if sys.dt ~= params.dt
                id = 'CCM:MustBeEqual';
                error(id,'The dynamics and the parameters have different sampling intervals defined. Please set them equal.')
            end
            
            % Check that number of disturbance states matches
            if sys.nw ~= params.nw
                id = 'CCM:MustBeEqual';
                error(id,'The dynamics and the parameters have different numbers of disturbance states defined. Please set them equal.')
            end
            
            % Check that setting for parametric uncertainty matches
            if sys.param_uncertainty ~= params.param_uncertainty
                id = 'CCM:MustBeEqual';
                error(id,'The dynamics and the parameters have a differently set parametric uncertainty. Please set them equal.')
            end
            % Enable all warnings and pause
            warning('on');
            pause('on')
            
            % Define system and parameters
            obj.sys = sys;
            obj.params = params;
            obj.use_sos = options.use_sos;
            obj.adaptive = options.adaptive;
                        
            % Define system dimension and horizon
            obj.nx = sys.nx;
            obj.nu = sys.nu;
            obj.np = sys.np;
            obj.nw = sys.nw;
            obj.N = params.N;

            % Define parametrization of CCM
            obj.W_idxs = idx2logical(options.W_idxs,obj.nx);
            obj.state2param = idx2mat(options.W_idxs,obj.nx);
            obj.monomial_degree = options.monomial_degree;
            
            % Define recursive feasibility and delta at time zero
            obj.terminal_constraint = options.terminal_constraint;
            obj.exact_initialization = options.exact_initialization;
            
            % Set number of grid points
            obj.n_tot = options.n_tot;
            obj.grid_pattern = options.grid_pattern;
            
            % Set booleans to enforce recomputation
            options.recompute_rccm = options.recompute.rccm || options.recompute_rccm;
            recompute = structfun(@(X) X || options.recompute_constants, options.recompute,'UniformOutput',false);
            recompute.rccm = options.recompute_rccm;
            
            % Define rccm variables and properties
            props = ["params.F_x","params.b_x", ...
                     "params.F_u","params.b_u", ...
                     "params.theta_v", "params.w_max", ...
                     "rho_des","W_idxs","monomial_degree","use_sos"];
            rccm_vars = ["rho","Y_fcn","Y_coef","W_fcn","W_coef","dW_dx","dW_dt_fcn","dW_dt_approx_fcn"];
            
            % Search optimal contraction rate rho_opt
            n_rho = numel(options.rho);
            tightenings = NaN(numel([obj.params.b_x; obj.params.b_u]),n_rho);
            for i = 1:n_rho
                % Set desired contraction rate
                obj.rho_des = options.rho(i);
                fprintf('Compute ccm with rho = %.3f\n',obj.rho_des);

                % Synthesise CCM
                [obj, success] = obj.loadVars('rccm',rccm_vars,props);
                if ~success || recompute.rccm
                    fprintf(2,'Computing rccm, this might take a while...\n')
                    [obj, success] = obj.synthesiseCCM(obj.rho_des);
                end
                
                % Compute constants
                if success
                    [obj, success] = obj.computeConstants(recompute);
                    fprintf('\n')

                    if success
                        obj.saveVars('rccm',rccm_vars,props);
                    end
                end
                
                % Compute asymptotic tightenings
                if success
                    tightenings(:,i) = [obj.c_x*obj.delta_inf;
                                        obj.c_u*obj.delta_inf];
                end
                
                % Notify user about invalid ccm
                if ~success
                    fprintf(2,'Warning: No valid ccm found for rho = %.3f\n',obj.rho_des);

                    if i == n_rho && all(isnan(tightenings),'all')
                        % Error that no valid ccm was found
                        id = 'CCM:FailedToConstructCCM';
                        error(id,'Failed to find a valid ccm.')
                    end
                end
            end
            % Load ccm with optimal contraction rate rho
            if n_rho > 1
                obj.rho_des = obj.getOptimalContractionRate(options.rho,tightenings);
                fprintf('Success, optimal contraction rate rho = %.3f found!\n',obj.rho_des)
                [obj, success] = obj.loadVars('rccm',rccm_vars,props);
                if success
                    [obj, success] = obj.computeConstants();
                end
                if ~success
                    error('Something went wrong. Optimal ccm could not be loaded.')
                end
            end
            
            % Setup geodesics
            obj.geod = setupGeodesic(obj.W_fcn,obj.dW_dx,obj.nx);
            
            % Define variable sizes
            obj.n_v = obj.N*obj.nu;
            obj.n_z = (obj.N+1)*obj.nx;
            obj.n_delta = obj.N+1;
            obj.n_w_bar = obj.N*4; % 4 because of intermediate rk4 points
            obj.n_geod = (obj.geod.D+1)*obj.nx;
            
            % Setup RAMPC solver
            obj = obj.setupRAMPC(obj.terminal_constraint, ...
                                 obj.exact_initialization, ...
                                 options.solver_options);

            % Initialize solver (using a nominal MPC)
            obj.y0 = obj.initializeSolver(obj.terminal_constraint, ...
                                          options.solver_options);
        end

        function [v, s, y0] = solve(obj,x0,theta0,options)
            arguments
                obj (1,1) CCM
                x0 (:,:) {mustBeVector,mustBeNumeric}
                theta0 (:,2) {mustBeNumeric} = obj.params.theta_v
                options.y0 (:,:) {mustBeNumeric,mustBeVector} = obj.y0
            end
            % Concatenate parameters
            par0 = [x0(:);reshape(theta0',[],1)];
            
            % Solve optimization
            tic;
            res = obj.solver('p',par0,'x0',options.y0(:),'lbx',obj.ylb,'ubx',obj.yub,'lbg',obj.clb,'ubg',obj.cub);
            s.solvetime = toc;
            
            % Get solver stats
            s.stats = obj.solver.stats();
            s.stats.feasible = s.stats.success;
            s.stats = rmfield(s.stats,'success');
            if ~s.stats.feasible
                warning(s.stats.return_status);
            end
            
            % Extract the solution
            s.y_sol = full(res.x);
            idx = obj.n_geod;
            s.v_sol = reshape(s.y_sol(idx + (1:obj.n_v)), [obj.nu,obj.N]);   idx = idx+obj.n_v;
            s.z_sol = reshape(s.y_sol(idx + (1:obj.n_z)), [obj.nx,obj.N+1]); idx = idx+obj.n_z;
            s.delta_sol = s.y_sol(idx + (1:obj.n_delta));                    idx = idx+obj.n_delta;
            if obj.adaptive
                s.theta_bar_sol = s.y_sol(idx+obj.n_w_bar+1);                idx = idx+obj.n_w_bar;
                s.delta_bar_sol = s.y_sol(idx+1);
            else
                s.theta_bar_sol = zeros(obj.np,1);
                s.delta_bar_sol = min(obj.delta_bar_x,obj.delta_bar_u);
            end
            v = s.v_sol(:,1);

            % Extract active constraints
            lam_g = full(res.lam_g);
            s.active_con_x = reshape(lam_g(obj.xineq_idxs), [numel(obj.params.b_x),obj.N+1]);
            s.active_con_u = reshape(lam_g(obj.uineq_idxs), [numel(obj.params.b_u),obj.N]);
            
            % Shift solver solution
            y0 = obj.shiftSolution(s.y_sol,s.delta_bar_sol,s.theta_bar_sol);
            
            % Extract cost
            s.cost = full(res.f);
            s.cost_regularizer = obj.params.regularizer*(s.delta_sol'*s.delta_sol);
            
            % Propagate tube scaling with exact equality
            s.delta_tight = obj.deltaProp(s.z_sol,s.v_sol,s.delta_sol,s.theta_bar_sol,theta0);
            
            % Extract tight tightenings
            s.tightenings_x = s.delta_tight'.*obj.c_x;
            s.tightenings_u = s.delta_tight(1:end-1)'.*obj.c_u;
            
            % Compute tubes
            s.tubes_x = getTubes(obj.params.F_x,s.tightenings_x);
            s.tubes_u = getTubes(obj.params.F_u,s.tightenings_u);
            
            % Define time vector
            s.t = 0:obj.params.dt:obj.params.T;
            
            % Save system inputs, offline computations and parameters
            s.inputs = struct();
            s.inputs.rho = obj.rho_des;
            s.inputs.W_idxs = obj.W_idxs;
            s.inputs.use_sos = obj.use_sos;
            s.inputs.adaptive = obj.adaptive;
            s.inputs.monomial_degree = obj.monomial_degree;
            s.inputs.terminal_invariance = obj.terminal_invariance;
            s.inputs.terminal_constraint = obj.terminal_constraint;
            s.inputs.exact_initialization = obj.exact_initialization;
            s.offline = struct();
            s.offline.rho = obj.rho;
            s.offline.Y_fcn = obj.Y_fcn;
            s.offline.W_fcn = obj.W_fcn;
            s.offline.c_x = obj.c_x;
            s.offline.c_u = obj.c_u;
            s.offline.c_obs = obj.c_obs;
            s.offline.L_D = obj.L_D;
            s.offline.L_G = obj.L_G;
            s.offline.M_under = obj.M_under;
            s.offline.delta_bar_x = obj.delta_bar_x;
            s.offline.delta_bar_u = obj.delta_bar_u;
            s.offline.delta_over = obj.delta_over;
            s.offline.delta_inf = obj.delta_inf;
            s.params = struct(obj.params);
        end
    end
    
    methods (Access=public)
        [obj, success] = synthesiseCCM(obj,rho)
        [obj, success] = solveSOS(obj,rho,v_W_fun,dv_W_dx_fun,W_coef,Y_coef)
        [obj, success] = solveLMIs(obj,rho,v_W_fun,dv_W_dx_fun,W_coef,Y_coef)
        rho = checkCCM(obj)
        [obj, success] = computeConstants(obj,recompute)
        obj = compute_c_x(obj,M_chol_fcn,hs_x_fcn);
        obj = compute_c_u(obj,M_chol_fcn,hs_u_fcn);
        obj = compute_c_obs(obj,x_sym,M_chol_fcn,hs_obs_fcn);
        obj = compute_L_D(obj,x_sym,d_sym,M_chol_fcn,Es_fcn);
        obj = compute_L_G(obj,x_sym,u_sym,theta_sym,M_chol_fcn,Gs_fcn);
        [d_delta, delta_inf] = checkConstants(obj,z_s0,v_s0,delta_bar_0,M_chol_fcn);
        obj = computeMunder(obj)
        rho = getOptimalContractionRate(OBJ,rho,tightenings);
        delta_tight = deltaProp(obj,z,v,delta,theta_bar,theta_v)
        ret = w_tilde(obj,x,u,theta_bar,theta,d,W_chol_x)
        y0 = shiftSolution(obj,y_sol,delta_bar_sol,theta_bar_opt)
        y0 = initializeSolver(obj,terminal_constraint,solver_options)
        u = kappa(obj,x,z,v)

        [obj, success] = loadVars(obj,name,vars,props)
        saveVars(obj,name,vars,props)
        s  =  property2struct(obj,props,s)
        obj = struct2property(obj,props,s)
        file_name = getFileName(obj)
    end
    
    methods (Access=private)
        TF = idx2mat(idxs,n)
        TF = idx2logical(idxs,n)
        TF = structcmp(s1,s2,exactMatch)
        TF = findDependencies(fcn,varargin)
        unique_field_name = getUniqueFieldName(s,field_name)
        field_name = findMatchingField(s1,s2,identifier,exactMatch)
        n_grid = getNgrid(obj,n_vars,grid_pattern,options)
        samples = getSamples(obj,F,b,n_grid,keep)
        comb = getCombinations(varargin)
        vert = poly2vertex(F,b)
        box = poly2box(obj,F_x,b_x,F_u,b_u,F_d,b_d,theta_v,options)
        [f_sym, sym_vars] = sdpvar2syms(varargin)
        poolobj = startParpool(profile)
        progbar(progress,options)
        
        [x,w] = clencurt(N)
        tubes = getTubes(F,tightenings)
        [T, T_dot] = computeCheby(N,D,t)
        E = RiemannEnergy(c,nx,D,N,T,T_dot,w,W_fcn)
        g = energyGradient(c,nx,D,N,T,T_dot,w,W_fcn,dW_fcn) 
        geodesic = setupGeodesic(W_fcn,dW_dx,nx)
        [gamma, gamma_s] = computeGeodesic(x,z,initial_geodesic)
        Vd = VdeltaSquared(N,gamma,gamma_s,w,W_fcn)
    end
end