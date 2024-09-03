function geodesic = setupGeodesic(W_fcn,dW_dx,nx)
    % Taken from Zhao https://github.com/boranzhao/robust_ccm_tube/blob/master/pvtol/sim/set_opt_prob_for_geodesic_computation.m
    
    % Problem setting for geodesic computation
    D = 2;      % Degree of the polynomial
    N = D+6;    % Stopping index for the CGL (Chebyshev-Gauss-Lobatto) nodes: #notes N+1
    
    % Optimization variables: chebyshev coefficients for geodesics
    % [c_10, c_11,..., c_1D, ..., c_n0, c_n1,..., c_nD]
    
    % Get chebyshev pseudospectral numerics (used for computing the integral)
    [s,w_cheby] = clencurt(N); % t is 1 by N+1, with values between 0 and 1.
    
    % Compute Cheby basis at all points
    [T, Tdot] = computeCheby(N,D,s); % Both T and T_dot are D+1 by N+1
    
    % Define equality constraints
    Aeq = [kron(eye(nx),T(:,1)'); kron(eye(nx),ones(1,D+1))];
    Aeq = sparse(Aeq);
    beq = zeros(2*nx,1);c0 = zeros(nx*(D+1),1);
    
    % Define cost and gradient
    ndec = nx*(D+1); %#ok
    costf = @(c) RiemannEnergy(c,nx,D,N,T,Tdot,w_cheby,W_fcn);
    grad  = @(c) energyGradient(c,nx,D,N,T,Tdot,w_cheby,W_fcn,dW_dx); 
    
    % Define lower and upper bounds on variables
    lb =-inf(nx,D+1);   % -20*ones(nx,D+1);
    ub = inf(nx,D+1);   %  20*ones(nx,D+1);
    lb = lb';lb = lb(:);
    ub = ub';ub = ub(:);
    
    % Define gradient object for fmincon
    gradObj = @(c) deal(costf(c), grad(c));
    
    % Fmincon solver settings
    nlprob.x0 = c0;
    nlprob.solver = 'fmincon';
    nlprob.algorithm = 'interior-point';
    nlprob.Aineq = [];
    nlprob.bineq = [];
    nlprob.Aeq = Aeq;
    nlprob.beq = beq;
    nlprob.lb = lb;
    nlprob.ub = ub;
    nlprob.objective = gradObj;
    nlprob.options = optimoptions('fmincon', ...
                                  Display='off', ...
                                  HessianApproximation='lbfgs',...
                                  MaxIterations=500, ...
                                  SpecifyObjectiveGradient=true, ...
                                  CheckGradients=false,...
                                  OptimalityTolerance=1e-4, ...
                                  FunctionTolerance=1e-4, ...
                                  FunValCheck='on', ...
                                  StepTolerance=1e-8);
    
    % Save parameters
    geodesic.D = D;
    geodesic.N = N;  % geodesic.ndec = ndec;
    geodesic.T = T;
    geodesic.Tdot = Tdot;
    geodesic.lb = lb;
    geodesic.ub = ub;
    geodesic.Aeq = Aeq;
    geodesic.beq = beq;
    geodesic.costf = costf;
    geodesic.grad = grad;
    geodesic.w_cheby = w_cheby;
    geodesic.nlprob = nlprob;
end