function [gamma, gamma_s] = computeGeodesic(x,z,initial_geodesic)
    % Taken from Zhao https://github.com/boranzhao/robust_ccm_tube
    persistent t_pre beq_pre copt_pre Erem_pre
    
    % Get dimensions
    nx = size(x,1);
    D = initial_geodesic.D; 
    
    % Intialize persisten variables
    if isempty(t_pre) % || (t == 0 && t_pre ~=0)
        t_pre = -3;
        beq_pre = zeros(2*nx,1);
        copt_pre = zeros(nx*(D+1),1);
        Erem_pre = Inf;
    end
    
    % Define linear equality constraints
    beq = [z;x];
    
    % Initialize values
    c0 = zeros(nx*(D+1),1);
    geodesic = initial_geodesic;
    
    % Vectorized format to improve computational efficiency
    i = 1:nx;    
    c0((i-1)*(D+1)+1,1) = z;
    c0((i-1)*(D+1)+2,1) = x - z;
    
    if norm(beq-beq_pre)<1e-8 && ~isinf(Erem_pre)
        copt = copt_pre;
        Erem = Erem_pre;
    else        
        % Solve optimization problem using fmincon
        geodesic.nlprob.x0 = c0;
        geodesic.nlprob.beq = beq;
        [copt,Erem,exitflag,~] = fmincon(geodesic.nlprob);
        if exitflag < 0
            warning('Geodesic optimization failed!');
        end
        
        % Update persistent variabels
        beq_pre = beq;
        copt_pre = copt;
        Erem_pre = Erem;
    end
    % Vectorized format (more computationally efficient)
    copt = transpose(reshape(copt,D+1,nx)); % the ith row corresponds to the ith element
    gamma = copt*geodesic.T;
    gamma_s = copt*geodesic.Tdot;
    
    % Compute control law
    % tic;
    % N = geodesic.N;
    % gamma = zeros(n,N+1);
    % gamma_s = zeros(n,N+1);  
    % for i = 1:n   
    %    gamma(i,:) = copt((i-1)*(D+1)+1:i*(D+1),:)'*T;       % gamma(i) is 1*(N+1); the ith elment of gamma on all the (N+1) nodes
    %    gamma_s(i,:) = copt((i-1)*(D+1)+1:i*(D+1),:)'*T_dot;
    % end  
    % toc;
    
    % Verify that the curve found is a geodesic (equation (11) Leung & Manchester)
    % error = 0;
    % for k=1:N+1
    %     error = error + (gamma_s(:,k)'*(controller.W_fcn(gamma(:,k))\gamma_s(:,k))-Erem)^2*geodesic.w_cheby(k);
    % end
    % error = sqrt(error)/Erem
    % if error >= 1e-5
    %      warning('Error = %.3e, the curve optimized is probably not a geodesic!\n',error);
    % end
end