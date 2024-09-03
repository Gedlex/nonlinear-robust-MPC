%% Simulate system
% clearvars;
% close all;
% clc
% 
% Load parameters
% params = PQRparams3();
% 
% % Create pre-stabilized dynamics
% sys = PlanarQuadrotor(params.dt,integrator=params.integrator,param_uncertainty=params.param_uncertainty);

% Initialize trajectory
N = 10000;
x = zeros(sys.nx,N+1);
x(:,1) = [-15,-6,0,0,0,0];
u_eq = ones(2,1)*sys.g*sys.m_nom/2;

% Simulate trajectory with u=u_eq
for i=1:N
    x(:,i+1) = sys.ddyn(x(:,i),u_eq);
end

disp(x(:,end))