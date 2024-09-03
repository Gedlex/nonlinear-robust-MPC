%% Search for Rho
% Load parameters
n_tot = 10;
grid_pattern = [1,1,1,1];
use_sos = false;
params = CPparams1();

% Create planar quadrotor
sys = CaptureStabilization(params.dt, ...
                           integrator=params.integrator, ...
                           approximate=use_sos, ...
                           param_uncertainty=params.param_uncertainty);

% Sample rho at origin
x = zeros(sys.nx,1); %x(3) = pi/4;
u = ones(sys.nu,1)*sys.g*sys.m_nom/2;
d = zeros(sys.nw,1);
theta = zeros(sys.np,1);
     
% Linearize dynamics
A = sys.A_diff(x,u,theta,d);
B = sys.B_diff(x,u,theta,d);

% Check controllability
Co = ctrb(A,B);
if length(A) == rank(Co)
    disp("The system is fully controllable.")
else
    error("The sytem has %d uncontrollable states.",length(A) - rank(Co))
end

% Compute LQR-gain
[~,~,P] = lqr(A,B,params.Q_cost,params.R_cost);

% Set rho based on slowest eigenvalue
rho = max(abs(P));

fprintf('Rho ~=: %.4f\n', rho);