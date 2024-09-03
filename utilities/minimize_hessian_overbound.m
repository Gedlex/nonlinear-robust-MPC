%% Minimize Hessian Overbound
clearvars;
close all;
clc;

% Load parameters
params = PQRparams2();
% params.integrator = "rkdp";
% params.dt = 0.1;

% Load system
sys = PlanarQuadrotor(params.dt,integrator=params.integrator,param_uncertainty=params.param_uncertainty);

%% Compute Hessian
import casadi.*
x_sym = SX.sym('x',sys.nx);
u_sym = SX.sym('u',sys.nu);
K_sym = SX.sym('K',[sys.nu,sys.nx]);
sym_vars = [x_sym; u_sym];

% Create stabilized dynamics
f = @(x,u) sys.f_fcn(x) + sys.B_fcn(x)*(K_sym*x + u);

% Compute multidimensional hessian
H_fcn = Function('H',{x_sym,u_sym,K_sym},{jacobian(jacobian(sys.ddyn(x_sym,u_sym,dynamics=f), sym_vars), sym_vars)});

% Clear symbolic variables
clear('x_sym','u_sym','K_sym','f')

%% Create Optimization Variables
import casadi.*
n_var = sys.nx + sys.nu;
y = MX.sym('y', n_var);
p = MX.sym('K', sys.nu*sys.nx);

% Extract variables
x_sym = y(1:sys.nx);
u_sym = y(sys.nx + (1:sys.nu));
K_sym = reshape(p,sys.nu,sys.nx);

% Evaluate Hessian
H_sym = H_fcn(x_sym,u_sym,K_sym);

%% Define Constraints
% Define state and input constraints
g_state = params.F_x*x_sym - params.b_x;
g_input = params.F_u*(K_sym*x_sym + u_sym) - params.b_u;

% Define constraints
g_ineq = [g_state;
          g_input];

% Define lower and upper bound on constraints
n_ineq = length(g_ineq);
lbg =-inf(n_ineq,1);
ubg = zeros(n_ineq,1);

% Preallocate x0 and objective
obj = zeros(1,sys.nx);
x0 = zeros(sys.nx,n_var);

%% Solve Optimization
tic
import casadi.*
for i = 1:sys.nx
    fprintf('Compute mu_%d\n',i);

    % Extract Hessian of x_k
    H_k = H_sym(i:sys.nx:end,:);

    % Skip if Hessian is zero
    if H_k.is_zero
        continue
    end
    
    % Define objective
    objective = -0.5*sum(sum(abs(H_k)));
    
    % Setup NPSOL
    nlp = struct('x', y, 'p',p, 'f', objective, 'g', g_ineq);
    opts = struct('ipopt', struct('max_iter',1000, 'print_level',0), 'print_time',0);
    solver = nlpsol('solver', 'ipopt', nlp, opts);

    % Define K
    K = zeros(sys.nu, sys.nx);
    
    % Solve Problem
    res = solver('x0',x0(i,:)','p',K(:),'lbg',lbg,'ubg',ubg);

    % Get solver stats
    stats = solver.stats();
    if ~stats.success
        warning(stats.return_status);
    end
    
    % Update initial guess and save solution
    x0(i,:) = full(res.x(1:n_var))';
    obj(i) = full(res.f);
end
toc

%% Check Solution
rerun_flag = false;
mu = -obj;
mu_check = zeros(1,sys.nx);

% Check mu
for i = 1:sys.nx
    % Extract x and u
    x = x0(i,1:sys.nx);
    u = x0(i,sys.nx + (1:sys.nu));
    
    % Evaluate Hessian
    H_i = H_fcn(x,u,K);
        
    % Check mu
    mu_i = eye(1,sys.nx);
    for j = 1:sys.nx
        H_k = full(H_i(j:sys.nx:end,:));
        
        mu_i(j) = 0.5*sum(sum(abs(H_k)));
        
        % Check if objective is (locally) optimal
        if mu_i(j) - mu(j) > 1E-5
            x0(j,:) = x0(i,:);
            rerun_flag = true;
        end
    end
    mu_check = max(mu_check,mu_i);
end