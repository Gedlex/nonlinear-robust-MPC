% -----------------------------------------------------------------------------
% File: SLS.m
% Author: Antoine Leeman (aleeman@ethz.ch)
% Date: 09th September 2023
% License: MIT
% Reference:
% Leeman AP., Sieber J., Bennani S., Zeilinger MN. “Robust optimal control 
% for nonlinear systems with parametric uncertainties via system level synthesis.” 
% Link: https://arxiv.org/abs/2304.00752
% -----------------------------------------------------------------------------

classdef SLS
    properties (GetAccess=public,SetAccess=private)
        sys (1,1) {mustBeA(sys,["NonlinearSystem","double"])} = NaN        % System dynamics
        params (1,1) {mustBeA(params,["Parameters","double"])} = NaN       % Parameters
        
        K (:,:) {mustBeNumeric}                                            % Pre-stabilizing feedback gain
        gamma_max (1,1) {mustBeNumeric}                                    % Performance index
        robust_performace_guarantees (1,1) logical                         % Boolean to activate rpgs

        N (1,1) {mustBeInteger,mustBeNonnegative}                          % Horizon (prediction steps)
        nx (1,1) {mustBeInteger,mustBeNonnegative}                         % Number of states
        nu (1,1) {mustBeInteger,mustBeNonnegative}                         % Number of inputs
        np (1,1) {mustBeInteger,mustBeNonnegative}                         % Number of parametric uncertainties
        nw (1,1) {mustBeInteger,mustBeNonnegative}                         % Number of disturbance states
                
        E (:,:) {mustBeNumeric}                                            % Discretized disturbance matrix
        mu (:,:) {mustBeNumeric,mustBeVector(mu,"allow-all-empties")}      % Overbound of the Lagrange error bound
        
        solver (1,1) {mustBeA(solver,["casadi.Function","double"])} = NaN  % NLP optimizer
    end

    properties (Access=public)
        y0 (:,:) {mustBeNumeric,mustBeVector(y0,"allow-all-empties")}      % Initial value of solver
    end

    properties (GetAccess=public,SetAccess=private,Hidden)
        n_tot_max (1,1) {mustBeInteger,mustBeNonnegative}                  % Number of maximally allowed grid points
        grid_pattern (4,1) {mustBeNonnegative}                             % Number of grid points relative to each other 
        n_eq (1,1) {mustBeInteger,mustBeNonnegative}                       % Number of equalities
        n_ineq (1,1) {mustBeInteger,mustBeNonnegative}                     % Number of inequalities
        n_y_nom (1,1) {mustBeInteger,mustBeNonnegative}                    % Number of nominal variables
        n_y_contr (1,1) {mustBeInteger,mustBeNonnegative}                  % Number of contraction variables
        n_y_tube (1,1) {mustBeInteger,mustBeNonnegative}                   % Number of tube variables 
        n_y_filter (1,1) {mustBeInteger,mustBeNonnegative}                 % Number of filter variables
        xineq_idxs (:,1) {mustBeInteger,mustBeNonnegative}                 % Indices of state inequality constraints
        uineq_idxs (:,1) {mustBeInteger,mustBeNonnegative}                 % Indices of input inequality constraints
        lbg (:,:) {mustBeNumeric,mustBeVector(lbg,"allow-all-empties")}    % Lower bounds on constraints
        ubg (:,:) {mustBeNumeric,mustBeVector(ubg,"allow-all-empties")}    % Upper bounds on constraints
    end

    methods
        function obj = SLS(sys,params,options)
            arguments
                sys (1,1) NonlinearSystem
                params (1,1) Parameters
                options.gamma_max (1,1) {mustBeNumeric} = 0.2
                options.recompute_E (1,1) logical = false
                options.recompute_mu (1,1) logical = false
                options.K (:,:) {mustBeNumeric} = zeros(sys.nu,sys.nx)
                options.robust_performace_guarantees (1,1) logical = true
                options.n_tot_max (1,1) {mustBeInteger,mustBePositive} = 1E6
                options.grid_pattern (4,1) {mustBeNonnegative} = [1,1,1,1]
                options.solver_options (1,1) struct = struct('ipopt', struct('max_iter',10000, 'print_level',5), 'print_time',0)
            end
            % Check for casadi
            if ~exist('casadi.MX', 'class')
                id = 'SLS:UndefinedClass';
                error(id,'CASADI not found. Add CASADI to the MATLAB search path.')
            end
            import casadi.*
            
            % Check for cprnd
            if ~exist('cprnd', 'file')
                id = 'SLS:UndefinedFunction';
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

            % Check prestabilizing K
            if (isprop(sys,'K') && any(options.K ~= sys.K,'all')) || (~isprop(sys,'K') && any(options.K,'all'))
                id = 'CCM:MustBeEqual';
                error(id,'The prestabilizing feedback for the dynamics must be the same as the one input to the SLS class. Please set them equal.')
            end
			% Enable all warnings and pause
            warning('on');
            pause('on')
            
            % Define system and parameters
            obj.sys = sys;
            obj.params = params;

            % Define system dimension and horizon
            obj.nx = sys.nx;
            obj.nu = sys.nu;
            obj.np = sys.np;
            obj.nw = sys.nw;
            obj.N = params.N;
            
            % Define prestabilizing controller and robust performance guarantees
            obj.K = options.K;
            obj.gamma_max = options.gamma_max;
            obj.robust_performace_guarantees =  options.robust_performace_guarantees;

            % Set number of grid points
            obj.n_tot_max = options.n_tot_max;
            obj.grid_pattern = options.grid_pattern;

            % Define properties to verify E and mu
            props = {'params.F_x','params.b_x', ...
                     'params.F_u','params.b_u', ...
                     'params.theta_v', ...
                     'params.dt','sys.integrator'};
            
            % Try to load mu
            [obj, success] = obj.loadVars('mu','mu',props);
            if ~success || options.recompute_mu
                fprintf(2,'Warning: Mu could not be loaded from file. Computing mu, this might take a while...\n')
                obj.startParpool();
                obj.mu = obj.compute_mu();
                obj.saveVars('mu','mu',props);
            end
            
            % Try to load E
            props = [props, 'params.w_max'];
            [obj, success] = obj.loadVars('E','E',props);
            if ~success || options.recompute_E
                fprintf(2,'Warning: E could not be loaded from file. Computing E, this might take a while...\n')
                obj.startParpool();
                obj.E = obj.compute_E();
                obj.saveVars('E','E',props);
            end
            
            % Initialize constraints
            obj.n_eq = 0;
            obj.n_ineq = 0;
            g_eq = [];
            g_ineq = [];
            var_slack = [];
            
            % Define optimization variables
            Z = MX.sym('State',obj.nx,obj.N+1);
            V = MX.sym('Input',obj.nu,obj.N+1);
            Phi_x = MX.sym('StateResponse',obj.nx,obj.nx,obj.N*(obj.N+1)/2);
            Phi_u = MX.sym('InputResponse',obj.nu,obj.nx,obj.N*(obj.N+1)/2);
            
            % Define parameter for initial state
            x0 = MX.sym('x0',obj.nx);
            
            % Define matrices for robust performance guarantees
            C = kron(eye(obj.N),eye(obj.nx));
            D = kron(eye(obj.N),eye(obj.nu));
            
            % Define filter and auxilary variables
            Sigma = cellmat(1,obj.N*(obj.N-1)/2,obj.nx,obj.nx);
            d = MX.sym('diag_Sigma',obj.N*obj.nx);
            tau = MX.sym('lin_tube',obj.N-1);
            
            % Get nominal variables
            [y_nom,obj.n_y_nom] = obj.getVariablesNominal(Z,V);
            
            % Constrain the nominal dynamics
            [n_dyn,g_dyn] = obj.getConstraintsDynamic(x0,Z,V);
            obj.n_eq = obj.n_eq + n_dyn;
            g_eq = [g_eq; g_dyn];
            
            % Nominal objective
            f = obj.getObjectiveNominal(Z,V);
            
            % Robust constraints satisfaction Eq. (35e)
            [n_cons,g_cons,slack_cons,n_slack_cons,xineq_idxs,uineq_idxs] = obj.getLinearConstraints(Z,V,Phi_x,Phi_u);
            obj.n_ineq = obj.n_ineq + n_cons;
            g_ineq = [g_ineq; g_cons];
            var_slack = [var_slack; slack_cons];
            
            % Filter constraints Eq. (35d)
            [n_filter,g_filter,n_slack_filter,slack_filter] = obj.getConstraintFilter(Z,V,Phi_x,Phi_u,d,tau);
            obj.n_ineq = obj.n_ineq + n_filter;
            g_ineq = [g_ineq; g_filter];
            var_slack = [var_slack; slack_filter];
            
            % System level parametrization of the controller Eq. (35b)
            [n_map,g_map] = obj.getNonlinearMapConstraints(Z,V,Phi_x,Phi_u,Sigma,d);
            obj.n_eq = obj.n_eq + n_map;
            g_eq = [g_eq; g_map];
            
            % Constraint on the tube size used for the linearization error overbound Eq. (35f)
            [n_ineq_tube,g_ineq_tube,n_slack_tube,slack_tube] = obj.getConstraintTube(tau,Phi_x,Phi_u);
            obj.n_ineq = obj.n_ineq + n_ineq_tube;
            g_ineq = [g_ineq; g_ineq_tube];
            var_slack = [var_slack; slack_tube];
            
            % Extract the variables
            [y_contr,obj.n_y_contr] = obj.getVariablesResponses(Phi_x,Phi_u);
            [y_tube,obj.n_y_tube] = obj.getVariablesTube(tau);
            [y_filter,obj.n_y_filter] = obj.getVariablesFilter_onlydiag(d);
            
            % Robust performance guarantees Eq.(42)
            n_var_slack_inf = 0;
            if obj.robust_performace_guarantees
                [inf_norm,n_ineq_inf,g_ineq_inf,var_slack_inf,n_var_slack_inf] = obj.mat_inf_norm([C*obj.v3_to_R(Phi_x);D*obj.v3_to_M(Phi_u)]);
                obj.n_ineq = obj.n_ineq + n_ineq_inf + 1;
                g_ineq = [g_ineq; g_ineq_inf; inf_norm-obj.gamma_max];
                var_slack = [var_slack; var_slack_inf; inf_norm];
            end
            
            % Shift idxs by equality constraints
            obj.xineq_idxs = xineq_idxs + obj.n_eq;
            obj.uineq_idxs = uineq_idxs + obj.n_eq;
            
            % Concatanate all the variables together
            y = [y_nom; y_contr; y_tube; y_filter; var_slack];
            
            % Get initial guess for solver variables
            [y_nom0,y_contr0,y_tube0,y_filter0,slack_filter0,slack_tube0] = obj.initializeSolver(y_nom,y_contr,y_tube,y_filter,slack_filter,slack_tube,g_map,g_filter,g_ineq_tube,options.solver_options);
            
            % Add a small regularization to the cost
            f = f + obj.params.regularizer*(y'*y);
            
            % Define lower and upper bound on constraints
            obj.lbg = [zeros(obj.n_eq,1);-inf(obj.n_ineq,1)];
            obj.ubg = zeros(obj.n_eq     +    obj.n_ineq, 1);
            
            % Initialize solver
            obj.y0 = [y_nom0;y_contr0;y_tube0;y_filter0;zeros(n_slack_cons,1);slack_filter0;slack_tube0;zeros(n_var_slack_inf,1)];
            
            % Initialize optimizer
            nlp = struct('x',y,'p',x0,'f',f,'g',[g_eq; g_ineq]);
            obj.solver = nlpsol('solver','ipopt',nlp,options.solver_options);
        end

        function [v, s, y0] = solve(obj,x0,options)
            arguments
                obj (1,1) SLS
                x0 (:,1) {mustBeNumeric}
                options.y0 (:,:) {mustBeNumeric,mustBeVector} = obj.y0
            end
            
            % Solve optimization
            tic;
            res = obj.solver('p',x0,'x0',options.y0(:),'lbg',obj.lbg,'ubg',obj.ubg);
            s.solvetime = toc;
            
            % Get solver stats
            s.stats = obj.solver.stats();
            s.stats.feasible = s.stats.success;
            s.stats = rmfield(s.stats,'success');
            if ~s.stats.feasible
                warning(s.stats.return_status);
            end
            
            % Return x as new initial guess
            y0 = full(res.x);
            
            % Extract the solution
            idx = 0;
            z_sol_v = full(res.x(idx+1:(obj.N+1)*obj.nx)); idx = idx+(obj.N+1)*obj.nx;
            v_sol_v = full(res.x(idx+1:obj.n_y_nom));      idx = obj.n_y_nom;
            Phi_sol_v = res.x(idx+(1:obj.n_y_contr));      idx = idx+obj.n_y_contr;
            tube_v = full(res.x(idx+(1:obj.n_y_tube)));    idx = idx+obj.n_y_tube;
            filter_v = full(res.x(idx+(1:obj.n_y_filter)));
            
            s.z_sol = obj.vecToNominalState(z_sol_v);
            s.v_sol = obj.K*s.z_sol + obj.vecToNominalInput(v_sol_v);
            [Phi_x_sol,Phi_u_sol] = obj.vecToResponse(Phi_sol_v);
            [Sigma_sol_v,d_sol] = obj.vecToSigma(filter_v);
            v = s.v_sol(:,1);
            
            % Transform Phi_u and Phi_x into Phi_u_tilde to account for feedback K
            Phi_u_tilde = cellfun(@(Phi_u,Phi_x) Phi_u + obj.K*Phi_x, Phi_u_sol, Phi_x_sol, 'UniformOutput', false);

            % Extract active constraints
            lam_g = full(res.lam_g);
            s.active_con_x = reshape(lam_g(obj.xineq_idxs), [numel(obj.params.b_x),obj.N+1]);
            s.active_con_u = reshape(lam_g(obj.uineq_idxs), [numel(obj.params.b_u),obj.N+1]);
            
            % Extract tubes
            s.R_sol = full(obj.v3_to_R(Phi_x_sol));
            s.M_sol = full(obj.v3_to_M(Phi_u_tilde));
            
            % Extract filter
            s.Sigma_sol = full(obj.v3_to_Sigma(Sigma_sol_v,d_sol));
            
            % Extract cost
            s.cost = full(res.f);
            s.cost_regularizer = obj.params.regularizer*(y0'*y0);
            
            % Compute tubes
            s.tubes_x = obj.getTubes(obj.nx,s.R_sol);
            s.tubes_u = obj.getTubes(obj.nu,s.M_sol);

            % Compute errors
            [s.B1,s.B2,s.B3] = obj.getErrors(d_sol,tube_v);
            
            % Define time vector
            s.t = 0:obj.params.dt:obj.params.T;

            % Save system inputs, offline computations and parameters
            s.inputs = struct();
            s.offline = struct();
            s.offline.E = obj.E;
            s.offline.mu = obj.mu;
            s.inputs.K = obj.K;
            s.inputs.gamma_max = obj.gamma_max;
            s.inputs.robust_performace_guarantees = obj.robust_performace_guarantees;
            s.params = struct(obj.params);
        end

        function [z_sol,v_sol,stats] = computeNominalMPC(obj,options)
            arguments
                obj (1,1) SLS
                options.solver_options (1,1) struct = struct('ipopt', struct('max_iter',1000, 'print_level',5), 'print_time',0)
            end
            import casadi.*
            Z = MX.sym('State',obj.nx,obj.N+1);
            V = MX.sym('Input',obj.nu,obj.N+1);
                       
            % Nominal objective
            f = obj.getObjectiveNominal(Z,V);
            
            % Constrain the nominal dynamics
            [n_eq,g_eq] = obj.getConstraintsDynamic(obj.params.x0,Z,V);%#ok
            
            % Nominal constraints (without robustness)
            [n_ineq, g_ineq, var_slack, ~] = obj.getLinearConstraints(Z,V,cellmat(obj.N*(obj.N+1)/2,obj.nx,obj.nx),cellmat(obj.N*(obj.N+1)/2,obj.nx,obj.nu)); %#ok
            
            % Concatanate all the variables together
            [y_nom, n_y_nom] = obj.getVariablesNominal(Z,V);           %#ok
            y = [y_nom; var_slack];

            % Define lower and upper bound on constraints
            lbg = [zeros(n_eq,1);-inf(n_ineq,1)];                      %#ok
            ubg =  zeros(n_eq    +    n_ineq,1 );                      %#ok

            % Solve the NLP
            nlp = struct('x',y,'f',f,'g',[g_eq;g_ineq]);
            solver = nlpsol('solver', 'ipopt', nlp, options.solver_options);%#ok
            
            % Solve problem
            res = solver('lbg',lbg,'ubg',ubg);                         %#ok
            stats = solver.stats();                                    %#ok
            
            % Extract the solution
            idx = 0;
            z_sol_v = full(res.x(idx+1:(obj.N+1)*obj.nx)); idx = idx+(obj.N+1)*obj.nx;
            v_sol_v = full(res.x(idx+1:n_y_nom));                      %#ok
            z_sol = obj.vecToNominalState(z_sol_v);
            v_sol = obj.vecToNominalInput(v_sol_v);
        end

        function [y_nom0,y_contr0,y_tube0,y_filter0,slack_filter0,slack_tube0] = initializeSolver(obj,y_nom,y_contr,y_tube,y_filter,slack_filter,slack_tube,g_map,g_filter,g_ineq_tube,solver_options)
            import casadi.*
            % Concatanate all the variables together
            y = [y_contr; y_tube; y_filter; slack_filter; slack_tube];
            
            % Define equality and inequality constraints
            g_eq = g_map;
            g_ineq = [g_filter; g_ineq_tube];

            % Get number of constraints
            n_eq = length(g_eq);                                       %#ok
            n_ineq = length(g_ineq);                                   %#ok
            n_slack_filter = length(slack_filter);
            
            % Compute objective
            f = (y'*y);
            
            % Define lower and upper bound on constraints
            lbg = [zeros(n_eq,1);-inf(n_ineq,1)];                      %#ok
            ubg =  zeros(n_eq    +    n_ineq, 1);                      %#ok
            
            % Initialize solver with all zeros
            y0 = zeros(length(y),1);                                   %#ok
            
            % Initialize optimizer
            nlp = struct('x', y, 'p', y_nom, 'f', f, 'g', [g_eq;g_ineq]);
            solver = nlpsol('solver', 'ipopt', nlp, solver_options);   %#ok
            
            % Solve nominal MPC for initial guess
            [z_sol,v_sol,stats] = obj.computeNominalMPC(solver_options=solver_options);
            
            % Get nominal variables
            [y_nom0,~] = obj.getVariablesNominal(z_sol,v_sol);
            if ~stats.success
                y_nom0 = 0*y_nom0;
            end
            
            % Solve problem
            res = solver('p',y_nom0,'x0',y0,'lbg',lbg,'ubg',ubg);      %#ok
            stats = solver.stats();                                    %#ok
            if stats.success
                y0 = res.x;                                            %#ok
            end

            % Extract solution
            idx = 0;
            y_contr0 = full(y0(idx+(1:obj.n_y_contr)));       idx = idx+obj.n_y_contr; %#ok
            y_tube0 = full(y0(idx+(1:obj.n_y_tube)));         idx = idx+obj.n_y_tube;  %#ok
            y_filter0 = full(y0(idx+(1:obj.n_y_filter)));     idx = idx+obj.n_y_filter;%#ok
            slack_filter0 = full(y0(idx+(1:n_slack_filter))); idx = idx+n_slack_filter;%#ok
            slack_tube0 = full(y0(idx+1:end));                                         %#ok
        end
        
        function R = v3_to_R(obj,v3_Phi_x)
            R = [];
            n=1;
            for k = 1:obj.N
                R_= [];
                for j=1:k
                    R_ = [R_, v3_Phi_x{n}];
                    n = n+1;
                end
                for j=k+1:obj.N
                    R_ = [R_, sparse(obj.nx,obj.nx)];
                end
                R = [R;R_];
            end
            
        end
        function M = v3_to_M(obj,v3_Phi_u)
            M = [];
            n=1;
            for k = 1:obj.N
                M_ = [];
                for j=1:k
                    M_ = [M_, v3_Phi_u{n}];
                    n = n+1;
                end
                for j=k+1:obj.N
                    M_ = [M_, sparse(obj.nu,obj.nx)];
                end
                M = [M;M_];
            end
        end
        
        function R = v3_to_Sigma(obj,v3_Sigma,d)
            R = [];
            n=1;
            nx = obj.nx;
            for k = 1:obj.N
                R_= [];
                for j=1:k-1
                    R_ = [R_, v3_Sigma{n}]; %%change here
                    n = n+1;
                end
                d_k = d((k-1)*nx+1:k*nx);
                R_ = [R_, diag(d_k)];
                for j=k+1:obj.N
                    R_ = [R_, sparse(nx,nx)];
                end
                R = [R;R_];
            end
            
        end
        
        
        function L = Phi_line_k(obj,k,v3_Phi)
            %k is the number of the line, starting from 1
            L = [];
            last = k*(k+1)/2;
            first = k*(k-1)/2+1;
            for j=first:last
                L = [L, v3_Phi{j}];
            end
        end
        
        function AB = buildNonlinearMap(obj,Z,V)
            nx = obj.nx;
            nu = obj.nu;
            N = obj.N;
            
            AB = sparse(nx, N*(nx+nu));
            AB(1:nx,1:nx) = eye(nx);
            
            for i =1:N-1
                A = obj.sys.A(Z(:,i+1),V(:,i+1));
                B = obj.sys.B(Z(:,i+1),V(:,i+1));
                mat_1 = sparse(nx,(i-1)*nx);
                mat_2 = [-A, eye(nx)];%-A(Z(:,i))
                mat_3 =  sparse(nx,(N-2)*nx + (i-1)*(nu-nx));
                mat_4 = -B;
                mat_5 = sparse(nx, (N-i)*nu );
                AB_ = [mat_1, mat_2, mat_3,mat_4,mat_5];
                AB = [AB;AB_];
            end
        end
        
        function J = getObjectiveNominal(obj,Z,V)            
            J =0;
            N = obj.N;
            for i=1:N
                J = J+(obj.K*Z(:,i) + V(:,i))'*obj.params.R_cost*(obj.K*Z(:,i) + V(:,i)) + (Z(:,i)-obj.params.x_ref)'*obj.params.Q_cost*(Z(:,i)-obj.params.x_ref);
            end            
            J = J+ (Z(:,N+1)-obj.params.x_ref)'*obj.params.P_cost*(Z(:,N+1)-obj.params.x_ref);
            
        end
        
        function [n,g] = getConstraintsDynamic(obj,x0,Z,V)
            g = [Z(:,1) - x0];
            N = obj.N;            
            
            for i=1:N
                g = [g; Z(:,i+1) - obj.sys.ddyn(Z(:,i),V(:,i))];
            end
            n = (N+1)*obj.nx;
        end
        
        function [n,g] = getNonlinearMapConstraints(obj,Z,V,Phi_x,Phi_u,Sigma, d)
            N = obj.N;
            nx = obj.nx;
            nu = obj.nu;
            
            R = obj.v3_to_R(Phi_x);
            M = obj.v3_to_M(Phi_u);
            Sigma_with_diag = obj.v3_to_Sigma(Sigma,d);
            
            AB = obj.buildNonlinearMap(Z,V);
            
            g_full = [reshape(AB*[R;M] - Sigma_with_diag, [N^2*nx^2,1])];
            [n,g] = obj.reducedMap(g_full);
        end
        
        function [n,g] = reducedMap(obj,g)
            N = obj.N;
            nx = obj.nx;
            v_NZ = find(reshape(kron(tril(ones(N),0),ones(nx)), [N^2*nx^2,1])); % remove all trivial equality constraints
            g = g(v_NZ);
            n = N*(N+1)/2*nx^2;
        end
        
        function [one_norm, n,g_ineq, slack, n_var] = one_norm(obj, v)
            import casadi.*
            n_slack = length(v);
            slack = MX.sym('slack',n_slack); % 1-norm always computed on vertical vectors
            one_norm = MX.sym('slack',1);
            
            g_ineq = [ v - slack; -v - slack; sum(slack) - one_norm];
            
            n = 2*n_slack +1;
            n_var = n_slack +1;
        end
        
        function [ one_norm_line, n_ineq,g_ineq, var_slack,n_var_slack ] = vec_one_norm(obj, M)            
            var_slack = [];
            n_var_slack = 0;
            n_ineq = 0;
            g_ineq = [];
            dim = size(M);
            one_norm_line = [];
            for i = 1 : dim(1)
                [one_norm,n,g, var, n_var] = obj.one_norm(M(i,:)');
                one_norm_line = [one_norm_line; one_norm];
                n_ineq = n_ineq+n;
                g_ineq = [g_ineq;g];
                var_slack = [var_slack;var];
                n_var_slack  = n_var_slack + n_var;
            end
        end
        
        function  [ inf_norm, n_ineq,g_ineq, var_slack,n_var_slack ] = mat_inf_norm(obj, M)
            import casadi.*
            inf_norm = MX.sym('slack_tube',1);
            
            var_slack = [];
            n_var_slack = 1;
            n_ineq = 0;
            g_ineq = [];
            dim = size(M);
            for i = 1 : dim(1)
                [one_norm,n,g, var, n_var] = obj.one_norm(M(i,:)');
                n_ineq = n_ineq+n;
                g_ineq = [g_ineq;g];
                var_slack = [var_slack;var;one_norm];
                n_var_slack  = n_var_slack + n_var;
                
                g_ineq = [g_ineq;one_norm - inf_norm];
                n_ineq = n_ineq + 1;
            end
        end
        
        function [n_ineq,g_ineq, var_cons, n_var_cons, xineq_idxs, uineq_idxs] = getLinearConstraints(obj,Z,V,v3_Phi_x,v3_Phi_u)
            F_x = obj.params.F_x;
            b_x = obj.params.b_x;
            F_u = obj.params.F_u;
            b_u = obj.params.b_u;
            
            var_cons = [];
            n_var_cons = 0;
            n_ineq = numel(b_x);
            xineq_idxs = (1:n_ineq)';
            g_ineq = [F_x*Z(:,1) - b_x]; %#ok
            for k=1:obj.N
                Phi_k = Phi_line_k(obj, k, v3_Phi_x);
                CPhi_x = F_x * Phi_k;                                
                [one_norm, n,g, var, n_var_slack] = vec_one_norm(obj, CPhi_x);         
                g_ineq = [g_ineq;g];
                n_ineq = n_ineq+n;
                var_cons = [var_cons;var;one_norm];
                n_var_cons = n_var_cons + n_var_slack;
                xineq_idxs = [xineq_idxs; n_ineq + (1:numel(b_x))'];
                g_ineq = [g_ineq; F_x*Z(:,k+1) + one_norm - b_x];
                n_ineq = n_ineq + length(b_x);
            end
            
            uineq_idxs = n_ineq + (1:numel(b_u))';
            g_ineq = [g_ineq;  F_u*(V(:,1) + obj.K*Z(:,1)) - b_u];
            n_ineq = n_ineq + length(b_u);
            for k=1:obj.N
                Phi_x_k = Phi_line_k(obj, k, v3_Phi_x);
                Phi_u_k = Phi_line_k(obj, k, v3_Phi_u);
                CPhi_u = F_u * (Phi_u_k + obj.K*Phi_x_k);
                [one_norm, n,g, var, n_var_slack] = vec_one_norm(obj, CPhi_u);         
                g_ineq = [g_ineq;g];
                n_ineq = n_ineq+n;
                var_cons = [var_cons;var;one_norm];
                n_var_cons = n_var_cons + n_var_slack;

                uineq_idxs = [uineq_idxs; n_ineq + (1:numel(b_u))'];
                g_ineq = [g_ineq; F_u*(V(:,k+1) + obj.K*Z(:,k+1)) + one_norm - b_u];
                n_ineq = n_ineq + length(b_u);
            end
        end

        
        function [n_ineq,g_ineq,n_slack, slack] = getConstraintFilter(obj,Z,V,v3_Phi_x,v3_Phi_u,d,delta)
            nx = obj.nx;
            v_P = vecnorm(obj.E,1,2);
            g_ineq = [];
            n_ineq = 0;
            
            slack = [];
            n_slack =0;
               
            for i=1:size(obj.params.theta_v,2)
                k =1 ;
                d_k = d((k-1)*nx+1:k*nx);
                
                theta = obj.params.theta_v(:,i);
                [one_norm, n,g, var, n_var_slack] = vec_one_norm(obj, obj.sys.ddyn_theta(Z(:,1), V(:,1))*theta);
                g_ineq = [g_ineq; g; one_norm + v_P - d_k];
                n_ineq = n_ineq + n +nx;
                
                slack = [slack;var;one_norm];
                n_slack = n_slack + n_var_slack;
                for k = 2:obj.N
                    d_k = d((k-1)*nx+1:k*nx);               
                    Phix_k = obj.Phi_line_k(k-1, v3_Phi_x);
                    Phiu_k = obj.Phi_line_k(k-1, v3_Phi_u);

                    LHS = [obj.sys.A_theta(Z(:,k),V(:,k))*theta*Phix_k + obj.sys.B_theta(Z(:,k), V(:,k))*theta*Phiu_k,...
                           obj.sys.ddyn_theta(Z(:,k), V(:,k))*theta]; 
                    [one_norm, n,g, var, n_var_slack] = vec_one_norm(obj, LHS);
                    g_ineq = [g_ineq;g];
                    n_ineq = n_ineq+n;

                    slack = [slack;var;one_norm];                
                    n_slack = n_slack + n_var_slack;

                    g_ineq = [g_ineq; one_norm + v_P + delta(k-1)^2*obj.mu' - d_k];
                    n_ineq = n_ineq + nx;
                end
            end
        end
        
    
        function [n_ineq,g_ineq, n_var_norm, var_norm] = getConstraintTube(obj,delta, v3_Phi_x, v3_Phi_u)
            g_ineq = [];
            n_ineq = 0;
            n_var_norm = 0;
            var_norm = [];
            
            for k = 1:obj.N-1
                Phix_k = obj.Phi_line_k(k, v3_Phi_x);
                Phiu_k = obj.Phi_line_k(k, v3_Phi_u);
                Phi_k = [Phix_k;Phiu_k];
                
                [ inf_norm, n_ineq_norm,g_ineq_norm, var_slack_norm,n_var_slack_norm ] = obj.mat_inf_norm(Phi_k);
                
                n_var_norm = n_var_norm + n_var_slack_norm;
                var_norm = [var_norm;var_slack_norm;inf_norm];
                n_ineq = n_ineq + n_ineq_norm;
                g_ineq = [g_ineq; g_ineq_norm];
                
                g_ineq = [g_ineq; inf_norm - delta(k)];
                n_ineq = n_ineq+1;
            end
        end
       
        
        function [y,n] = getVariablesNominal(obj,Z,V)
            N = obj.N;
            nx = obj.nx;
            nu = obj.nu;
            y = [reshape(Z,[(N+1)*nx,1]); reshape(V,[(N+1)*nu,1])];
            n = (N+1)*(nx+nu);
        end
        
        function [y,n] = getVariablesTube(obj,delta)
            N = obj.N;
            y = reshape(delta,[N-1,1]);
            n = N-1;
        end
        
        function [y,n] = getVariablesResponses(obj,Phi_x,Phi_u)
            y = [];
            N = obj.N;
            nu = obj.nu;
            nx = obj.nx;
            for i = 1:(N+1)*N/2
                y = [y; reshape(Phi_u{i},[nu*nx,1])];
            end
            
            for i = 1:(N+1)*N/2
                y = [y;reshape(Phi_x{i},[nx^2,1])];
            end
            n = (N+1)*N/2 * nu*nx + (N+1)*N/2*nx^2;
        end
        
       
        function [y,n] = getVariablesFilter_onlydiag(obj,d)
            y = [];
            N = obj.N;
            nu = obj.nu;
            nx = obj.nx;
            y = [d];
            n = N*nx;
        end
         
        function Z = vecToNominalState(obj,y)
            nx = obj.nx;
            N = obj.N;
            Z = reshape(y,[nx,N+1]);
        end
        
        function V = vecToNominalInput(obj,y)
            nu = obj.nu;
            N = obj.N;
            V = reshape(y,[nu,N+1]);
        end
       
        function [Phi_x,Phi_u] = vecToResponse(obj,y)
            nu = obj.nu;
            N = obj.N;
            nx = obj.nx;
            k = 1;
            for n =1: (N+1)*N/2
                Phi_u{n} = full(reshape(y(k:k+nu*nx-1),[nu,nx]));
                k = k + nu*nx;
            end
            for n =1: (N+1)*N/2
                Phi_x{n} = full(reshape(y(k:k+nx*nx-1),[nx,nx]));
                k = k + nx^2;
            end
        end
        
        function [Sigma, d] = vecToSigma(obj,y)
           nu = obj.nu;
           N = obj.N;
           nx = obj.nx;
           k=1;
           d = y(k : k+nx*N-1);
           Sigma = cellmat(1,N*(N-1)/2,nx,nx);
        end

        function [tubes] = getTubes(obj,n,R)
            tubes = reshape(vecnorm(kron(eye(obj.N), eye(n)) * R,1,2), [n,obj.N]);
            tubes = [zeros(n,1), tubes];
        end

        function [B1,B2,B3] = getErrors(obj,d_sol,tube_v)
            B1 = abs(reshape(kron([0;tube_v].^2,obj.mu'),[obj.nx,obj.N]));                % contribution of lin error
            B2 = abs(reshape(kron(ones(obj.N,1),vecnorm(obj.E,1,2)),[obj.nx,obj.N]));     % contribution of add. noise
            B3 = abs(reshape(d_sol,[obj.nx,obj.N]));                                      % total tube
        end

        function [E_max] = compute_E(obj) % Estimation of the value of E via sampling
            % Determine variables used / needed in fcn
            fcn = @(x,u,theta,d) obj.sys.ddyn(x,u,theta=theta,d=d) - obj.sys.ddyn(x,u,theta=theta);
            vars_used = obj.findDependencies(fcn,obj.nx,obj.nu,obj.np,obj.nw);
            
            % Compute n_vars
            n_vars = cellfun(@sum, vars_used);
            
            % Compute number of grid points according to grid pattern
            n_grid = obj.getNgrid(n_vars,obj.grid_pattern,n_tot_max=obj.n_tot_max);
            fprintf(['The number of grid points used:       State %d\n', ...
                     '                                      Input %d\n', ...
                     '                                      Theta %d\n', ...
                     '                                Disturbance %d\n'], n_grid)
            
            % Get samples from constraint polyhedrons
            x_eval = obj.getSamples(obj.params.F_x,obj.params.b_x,n_grid(1),vars_used{1});
            u_eval = obj.getSamples(obj.params.F_u,obj.params.b_u,n_grid(2),vars_used{2});
            theta_eval = obj.getSamples(obj.params.F_theta,obj.params.b_theta,n_grid(3),vars_used{3});
            d_eval = obj.getSamples(obj.params.F_d,obj.params.b_d,n_grid(4),vars_used{4});
            n_eval = size(x_eval,2);
            
            % Initialize E_max
            E_max = zeros(obj.nx,1);
            
            % Create additional variables to prevent communication overhead (for parfor)
            ddyn_w = @(x,u,theta,d) obj.sys.ddyn(x,u,theta=theta,d=d);
            ddyn = @(x,u,theta,d) obj.sys.ddyn(x,u,theta=theta);
            K = obj.K;       %#ok
            F = obj.params.F_x;
            b = obj.params.b_x;
            
            % Create dataqueue for parfor progbar
            q = parallel.pool.DataQueue;
            afterEach(q, @(~) obj.progbar(increment=100/n_eval));
            
            % Initialize progress bar
            obj.progbar(0,prefix='Compute E: ')

            tic
            % Loop over state samples
            parfor (i = 1:n_eval)
                % Get current state
                x = x_eval(:,i);
                
                % Loop over input samples
                for j = 1:size(u_eval,2)
                    u = u_eval(:,j);
                    
                    % Compute shifted input
                    v = u - K*x;       %#ok
                    
                    % Loop over theta
                    for k = 1:size(theta_eval,2)
                        theta = theta_eval(:,k);

                        for l = 1:size(d_eval,2)
                            d = d_eval(:,l);
                            
                            % Evaluate nonlinear discretization with noise
                            x_nl = ddyn_w(x,v,theta,d);
                            
                            % Reject x_nl outside of the state constraints
                            if ~all(F*x_nl <= b)
                                continue
                            end
                            
                            % Compute difference of true and nominal discretized dynamics
                            E_tmp = x_nl - ddyn(x,v,theta);
                            E_max = max(E_max, abs(E_tmp));
                        end
                    end
                end
                % Update progress bar
                send(q,i);
            end
            E_max = diag(E_max);

            % Delete dataqueue
            delete(q)
            toc
        end

        function [mu_max] = compute_mu(obj) % Estimation of the value of mu via sampling
            % Determine variables used / needed in Hessian
            vars_used = obj.findDependencies(obj.sys.H,obj.nx,obj.nu,obj.np);
            
            % Compute n_vars
            n_vars = cellfun(@sum, vars_used);
            
            % Compute number of grid points according to grid pattern
            n_grid = obj.getNgrid(n_vars,obj.grid_pattern(1:3),n_tot_max=obj.n_tot_max);
            fprintf(['The number of grid points used:       State %d\n', ...
                     '                                      Input %d\n', ...
                     '                                      Theta %d\n'], n_grid)
            
            % Get samples from constraint polyhedrons
            x_eval = obj.getSamples(obj.params.F_x,obj.params.b_x,n_grid(1),vars_used{1});
            u_eval = obj.getSamples(obj.params.F_u,obj.params.b_u,n_grid(2),vars_used{2});
            theta_eval = obj.getSamples(obj.params.F_theta,obj.params.b_theta,n_grid(3),vars_used{3});
            n_eval = size(x_eval,2);
            
            % Initialize mu_max
            mu_max = -inf(1,obj.nx);
            
            % Create additional variables to prevent communication overhead (for parfor)
            H = @(x,u,theta) full(obj.sys.H(x,u,theta));
            K = obj.K;   %#ok
            nx = obj.nx; %#ok
            nu = obj.nu; %#ok
            
            % Create dataqueue for parfor progbar
            q = parallel.pool.DataQueue;
            afterEach(q, @(~) obj.progbar(increment=100/n_eval));
            
            % Initialize progress bar and dataqueue (used for parfor)
            obj.progbar(0,prefix='Compute mu: ')
            
            tic
            % Loop over state samples
            parfor (i = 1:n_eval)
                % Get current state
                x = x_eval(:,i);
                
                % Loop over input samples
                for j = 1:size(u_eval,2)
                    u = u_eval(:,j);

                    % Compute shifted input
                    v = u - K*x;       %#ok
                    
                    % Loop over theta
                    for k = 1:size(theta_eval,2)
                        theta = theta_eval(:,k);
                        
                        % Evaluate Hessian
                        H_k = permute(reshape(H(x,v,theta),[nx,nx+nu,nx+nu]),[3,2,1]); %#ok
                        
                        % Compute mu
                        d = size(H_k,3);
                        mu_k = eye(1,d);
                        for l = 1:d
                            mu_k(l) = 0.5*sum(sum(abs(H_k(:,:,l))));
                        end
                        
                        % Find maximum on hessian
                        mu_max = max(mu_max, mu_k);
                    end
                end
                % Update progress bar
                send(q,i);
            end
            toc
            % Delete dataqueue
            delete(q)
        end
    end

    methods % (Access=protected,Hidden)
        function samples = getSamples(obj,F,b,n_grid,keep)
            arguments
                obj (1,1) SLS
                F (:,:) {mustBeNumeric}
                b (:,:) {mustBeNumeric,mustBeVector}
                n_grid (1,1) {mustBeInteger,mustBePositive}
                keep (:,:) logical {mustBeVector} = true(size(F,2),1)
            end
            
            % Convert polytope to vertices
            vertices = obj.poly2vertex(F,b);
            n = size(F,2);
            
            % Check if conversion succeeded
            if ~isempty(vertices)
                % Initialize zero vector
                samples = repmat({0},n,1);
                
                % Grid over vertices in keep
                for k = 1:n
                    if keep(k)
                        samples{k} = linspace(vertices(k,1),vertices(k,2),n_grid);
                    end
                end
                
                % Permute vector to get all grid combinations
                samples = obj.getCombinations(samples{:});
            else
                % Inform user of possible inaccuracies
                warning(['Polytope can`t be separated and must be sampled at random. ' ...
                         'This might be less accurate than gridding. ' ...
                         'Check the result carefully or consider defining polytopes that can be gridded.'])
        
                % Fix random seed for reproducable results
                rng(0,'twister');
                
                % Get unconstrained states
                constr = any(F,1)';
                
                % Compute number of samples (to be equivalent to 2*gridding)
                n_sample = 2*n_grid^(sum(keep(:) & constr));
                
                % Initialize zero vector
                samples = zeros(n,n_sample);
                
                % Sample from constrained polyhedron only
                samples(constr,:) = cprnd(n_sample, F(:,constr), b)';
            end
        end

        function n_grid = getNgrid(obj,n_vars,grid_pattern,options)
            arguments
                obj (1,1) SLS
                n_vars (:,:) {mustBeVector,mustBeInteger}
                grid_pattern (:,:) {mustBeNumeric,mustBeNonnegative}
                options.n_tot_max (1,1) {mustBeInteger,mustBePositive} = 1E6
                options.n_grid (:,:) {mustBeVector,mustBeInteger,mustBeNonnegative} = zeros(numel(n_vars),1)
            end
            % Check inputs
            mustBeMember(numel(n_vars),numel(grid_pattern));
            mustBeMember(numel(n_vars),numel(options.n_grid));
            
            % Initialize n_grid
            n_grid = options.n_grid(:);
            
            % Normalize grid pattern
            grid_pattern = grid_pattern./grid_pattern(1);
            
            % Normalize n_max by already defined grid parameters
            idxs = (n_grid==0);
            if ~all(idxs)
                options.n_tot_max = options.n_tot_max / prod(n_grid(~idxs).^n_vars(~idxs));
            end
            
            % Compute number of normalized gridding points
            n = (options.n_tot_max / prod(grid_pattern(idxs).^n_vars(idxs)))^(1/sum(n_vars(idxs)));
            
            % Compute n_grid
            n_grid(idxs) = n.*grid_pattern(idxs);
            n_grid(n_vars==0) = 1;
            
            % Find all possible combinations
            n_grid_comb = obj.getCombinations([floor(n_grid), ceil(n_grid)]);
            
            % Compute all possible n_tot_max
            n_tot_comb = prod(n_grid_comb.^n_vars,1);
            
            % Find minimum
            [~,idx] = min(abs(n_tot_comb - options.n_tot_max));
            
            % Return n_grid being closest to n_tot_max
            n_grid = n_grid_comb(:,idx);
        end
        
        function saveVars(obj,name,vars,props)           
	        % Check if file already exists
            filepath = obj.getFileName();
            lock = [filepath,'.lock'];
            if isfile(filepath)
                % Wait for lock
                while isfile(lock)
                    pause(rand)
                end
                
                % Create lock
                fid = fopen(lock,'w');
                fclose(fid);
                
                % Load structured array
                data = load(filepath);
                
                % Delete lock
                delete(lock)

                % Convert properties to structured array
                props_struct = obj.property2struct(props);
                
                % Find field with prefix name
                field_name = obj.findMatchingField(data,props_struct,name);
                
                % Check if field was found
                props = [vars, props];
                if isempty(field_name)
                    % Create unique field name
                    field_name = obj.getUniqueFieldName(data,name);
                    
                    % Create new struct from properties
                    s.(field_name) = obj.property2struct(props);
                else
                    % Add properties to existing struct
                    s.(field_name) = obj.property2struct(props,data.(field_name)); 
                end
                
                % Wait for lock
                while isfile(lock)
                    pause(rand)
                end
                
                % Create lock
                fid = fopen(lock,'w');
                fclose(fid);
                
                % Save offline computations
                save(filepath,'-struct','s',field_name,'-append');
                
                % Delete lock
                delete(lock)
            else
                % Create lock
                fid = fopen(lock,'w');
                fclose(fid);
                
                % Save offline computations
                field_name = [name,'_0'];
                props = [vars, props];
                s.(field_name) = obj.property2struct(props);
                save(filepath,'-struct','s',field_name);
                
                % Delete lock
                delete(lock)
            end
        end

        function [obj, success] = loadVars(obj,name,vars,props)
            % Initialize return
            success = false;
            
            % Check if file exists
            filepath = obj.getFileName();
            lock = [filepath,'.lock'];
            if isfile(filepath)
                % Wait for lock
                while isfile(lock)
                    pause(rand)
                end
                
                % Create lock
                fid = fopen(lock,'w');
                fclose(fid);
                
                % Load structured array
                data = load(filepath);
                
                % Delete lock
                delete(lock)
                
                % Convert properties to structured array
                props = obj.property2struct(props);
                
                % Find field with prefix name that matches all properties in props exactly
                field_name = obj.findMatchingField(data,props,name,true);
                
                % Load vars from file
                if ~isempty(field_name)
                    try % to load variables from file
                        obj = obj.struct2property(vars,data.(field_name));
                        success = true;
                    catch
                        return
                    end
                end
            end
        end

        function field_name = findMatchingField(obj,s1,s2,prefix,exactMatch)
            % Get struct fieldnames and allocate return value
			fields = fieldnames(s1);
            field_name = '';

            % Boolean for non-existing fields
            if ~exist('exactMatch','var')
                exactMatch = false;
            end
            
            % Loop through struct and search for matching entry
            for i = 1:numel(fields)
                % Check only field names that start with the correct prefix
                if ~startsWith(fields{i},prefix)
                    continue
                end
				% Check if structs match
                if ~obj.structcmp(s2,s1.(fields{i}),exactMatch)
					continue
                end
				% Return field name found
				field_name = fields{i};
				break
            end
        end

        function tf = structcmp(obj,s1,s2,exactMatch)
            % Get struct fieldnames and allocate return value
            fields = fieldnames(s1);
            tf = false;
            
            % Compare structs against each other
            for i = 1:numel(fields)
                % Check if field exists
                if ~isfield(s2, fields{i})
                    if exactMatch
                        return
                    else
                        continue
                    end
                end
                % Recursively call structcmp to compare nested structs
                if isstruct(s1.(fields{i})) && isstruct(s2.(fields{i}))
                    if ~obj.structcmp(s1.(fields{i}),s2.(fields{i}),exactMatch)
                        return
                    end
                % Use isequal to compare other data
                elseif ~isequal(s1.(fields{i}), s2.(fields{i}))
                    return
                end
            end
            tf = true;
        end
        
        function s = property2struct(obj,props,s)
            % Convert input to character cell array
            props = cellstr(props);
            
            % Split dynamic properties into substrings for get- and setfield
            props = cellfun(@(s) split(s,'.'), props,'UniformOutput',false);
            
            % Allocate output struct
            if ~exist('s','var')
                s = struct();
            end
            
            % Loop through properties and append them to the struct
            for i = 1:numel(props)
                try % to set value
                    prop = getfield(obj, props{i}{:});
                    s = setfield(s, props{i}{:}, prop);
                catch
                    id = 'SLS:InvalidProperty';
                    error(id,'%s is not a valid property',props{i}{end})
                end
            end
        end

        function obj = struct2property(obj,props,s)
            % Convert input to character cell array
            mustBeText(props)
            props = cellstr(props);
            
            % Split dynamic properties into substrings for get- and setfield
            props = cellfun(@(s) split(s,'.'), props,'UniformOutput',false);
            
            % Loop through struct and assign properties from struct to class
            for i = 1:numel(props)
                try % to set value
                    prop = getfield(s, props{i}{:});
                    obj = subsasgn(obj, struct('type',repmat({'.'},numel(props{i}),1),'subs',props{i}), prop);
                catch
                    id = 'CCM:InvalidFieldorProperty';
                    error(id,'%s is not a valid field or property',strjoin(props{i},'.'))
                end
            end
        end
        
        function file_name = getFileName(obj)
            % Extract class directory
            [path, ~, ~] = fileparts(mfilename('fullpath'));
            
            % Return file name
            file_name = fullfile(path,[class(obj.sys),'.mat']);
        end
    end

    methods (Static)
        function unique_field_name = getUniqueFieldName(s,field_name)
            suffix = 0;
            % Test suffix until unique
            unique_field_name = sprintf('%s_%d',field_name,suffix);
            while isfield(s, unique_field_name)
                unique_field_name = sprintf('%s_%d',field_name,suffix);
                suffix = suffix + 1;
            end
        end

        function TF = findDependencies(fcn,varargin)
            import casadi.*
            % Check function
            n = numel(varargin);
            if isa(fcn,'function_handle') || isa(fcn,'casadi.Function')
                % Check input
                cellfun(@(X) mustBeNumeric(X), varargin);
                
                % Create symbolic variables
                sym_vars = cell(n,1);
                for i = 1:n
                    % Get size from input
                    inp = varargin{i};
                    if isscalar(inp)
                        dim = [inp,1];
                    elseif numel(inp) == 2
                        dim = inp(:);
                    else
                        dim = size(inp);
                    end
                    
                    % Create symbolic variable
                    name = sprintf('x%u',i);
                    sym_vars{i} = SX.sym(name,dim);
                end
                
                % Transform function to symbolic expression
                sym_fcn = fcn(sym_vars{:});
            else
                % Check input
                mustBeA(fcn,["casadi.SX","casadi.MX","casadi.DM"]);
                cellfun(@(X) mustBeA(X,["casadi.SX","casadi.MX","casadi.DM"]), varargin);
                
                % Pass inputs on
                sym_fcn = fcn;
                sym_vars = varargin;
            end
            
            % Determine symbolic variables used in sym_fcn
            TF = cell(n,1);
            for i = 1:n
                TF{i} = false(size(sym_vars{i}));
                for j = 1:numel(sym_vars{i})
                    TF{i}(j) = depends_on(sym_fcn, sym_vars{i}(j));
                end
            end
        end
        
        function comb = getCombinations(varargin)
            % Convert stacked vectors into cell arrays
            cellvec = cellfun(@(X) mat2cell(X, ones(1,size(X,1))), varargin, 'UniformOutput', false)';
            cellvec = vertcat(cellvec{:}); % flatten nested cells
            
            % Remove duplicates from vector
            cellvec = cellfun(@(X) unique(X), cellvec, 'UniformOutput', false);
            
            % Compute n-dimensional grid
            out = cell(size(cellvec,1),1);
            [out{:}] = ndgrid(cellvec{:});
            
            % Convert ndgrid into permuation matrix
            out = cellfun(@(X) X(:)', out, 'UniformOutput', false);
            comb = cell2mat(out);
        end
        
        function vert = poly2vertex(F,b)
            % Check boundedness for constrained states
            vert = [];
            if any((~any(F < 0) | ~any(F > 0)) & any(F))
                warning('Polytope %s is unbounded for constraint state.',inputname(1));
                return
            % Check separability
            elseif any(sum(F ~= 0, 2) > 1)
                warning('Polytope %s cannot be separated into simple 1D boxes.',inputname(1));
                return
            end
            
            % Compute 1D box constraints
            n = size(F,2);
            vert = zeros(n,2);
            for i = 1:n
                % Skip if state is unconstrained
                if ~any(F(:,i))
                    continue
                end
                
                % Find nonzero entries
                idx_p = (F(:,i) > 0);
                idx_n = (F(:,i) < 0);
                
                % Find infimum and supremum
                maximum = min(b(idx_p) ./ F(idx_p,i));
                minimum = max(b(idx_n) ./ F(idx_n,i));
                
                % Check that set is nonempty
                if minimum > maximum
                    error('Constraint on state %s is empyt set.',num2str(i));
                end
                
                % Add vertices to vert
                vert(i,:) = [minimum, maximum];
            end
        end

        function poolobj = startParpool(profile)
            arguments
                profile (1,1) string = parallel.defaultClusterProfile
            end
            poolobj = gcp("nocreate"); % If no pool, do not create new one.
            if isempty(poolobj)
                poolobj = parpool(profile);
            elseif ~strcmp(poolobj.Cluster.Profile,profile)
                delete(poolobj)
                poolobj = parpool(profile);
            end
        end

        function progbar(progress,options)
            arguments
                progress (1,1) {mustBeInRange(progress,0,100)} = 0
                options.increment (1,1) {mustBeNonnegative} = 0
                options.barLength (1,1) {mustBeInteger} = 0
                options.reset (1,1) logical = 0
                options.prefix (1,1) string = ""
            end
            persistent storedProgress prefix barLength
            
            % Initialize persistent variables
            if isempty(storedProgress) || options.reset
                storedProgress = 0;
            end
            if isempty(prefix) || options.reset
                prefix = 'Progress: ';
                barLength = 25;
            end
            
            % Update by increment
            if options.increment && ~progress
                progress = min(storedProgress + options.increment, 100);
            end
            
            % Generate reverse string
            strLength = 0;
            if 0 < progress && storedProgress <= progress
                strLength = strlength(prefix) + barLength + 8;
            end
            reverseStr = repmat(sprintf('\b'), 1, strLength);
            
            % Change persistent variables
            if options.prefix ~= ""
                prefix = options.prefix;
            end
            if options.barLength
                barLength = options.barLength;
            end
            
            % Build bar
            n_bars = round(progress*barLength/100);
            bar = sprintf('[%s%s]',repmat('#',1,n_bars),repmat('-',1,barLength-n_bars));
            
            % Print progress
            msg = sprintf('%s%s %3u%% ',prefix,bar,round(progress));
            fprintf([reverseStr,'%s'],msg);

            % Check for end
            if (100-progress < 1E-5) && (100-storedProgress > 1E-5)
                fprintf('\n');
            end

            % Update stored progress
            storedProgress = progress;
        end
    end
end

