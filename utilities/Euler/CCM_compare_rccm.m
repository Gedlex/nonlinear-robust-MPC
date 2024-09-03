%% Compare RCCMs
clearvars
clc

%% Compute RCCM from Paper
% Load parameters
[params, obs] = PQRparams_theta();

% Create planar quadrotor
sys = PlanarQuadrotor(params.dt, ...
                      integrator=params.integrator, ...
                      param_uncertainty=params.param_uncertainty);

% Create ccm class instance
ccm = CCM(sys,params,adaptive=false, ...
          exact_initialization=false, ...
          terminal_constraint=true, ...
          use_sos=true);

% Create table for constants
sz = [length(ccm.c_x) 8];
varTypes = repmat("double",1,8);
varNames = ["c_x","c_u","L_D","L_G","delta bar x","delta bar u","delta over","rho"];
T = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);

% Populate table
T{1:length(ccm.c_x),1} = ccm.c_x;
T{1:length(ccm.c_u),2} = ccm.c_u;
T{1:length(ccm.L_D),3} = ccm.L_D;
T{1:length(ccm.L_G),4} = ccm.L_G;
T{1:length(ccm.delta_bar_x),5} = ccm.delta_bar_x;
T{1:length(ccm.delta_bar_u),6} = ccm.delta_bar_u;
T{1:length(ccm.delta_over),7} = ccm.delta_over;
T{1:length(ccm.rho),8} = ccm.rho;

% Print constants
disp('Original Constants:')
disp(T)

%% Compute RCCM using LMIs
% Create struct defining sample numbers
n_tot = struct('rccm',15E3,'rccm_check',1E5,'c_x',1E5,'c_u',1E5,'c_obs',1E5,'L_D',1E5,'L_G',1E5,'w_bar',1E5,'delta_over',1E5,'M_under',5E3,'M_under_check',1E4);

n = [1E3,15E3,20E3,25E3,30E3,50E3];

for i = 1:numel(n)
    % Set number of LMI samples
    n_tot.rccm = n(i);
    
    % Create ccm class instance
    ccm = CCM(sys,params,adaptive=false, ...
              exact_initialization=false, ...
              terminal_constraint=false, ...
              recompute_rccm=true, ...
              use_sos=false, ...
              n_tot=n_tot);
    
    % Populate table
    T{:,:} = 0;
    T{1:length(ccm.c_x),1} = ccm.c_x;
    T{1:length(ccm.c_u),2} = ccm.c_u;
    T{1:length(ccm.L_D),3} = ccm.L_D;
    T{1:length(ccm.L_G),4} = ccm.L_G;
    T{1:length(ccm.delta_bar_x),5} = ccm.delta_bar_x;
    T{1:length(ccm.delta_bar_u),6} = ccm.delta_bar_u;
    T{1:length(ccm.delta_over),7} = ccm.delta_over;
    T{1:length(ccm.rho),8} = ccm.rho;
    
    % Print constants
    fprintf('RCCM with %.1E LMI-Samples:\n', n(i))
    disp(T)
end