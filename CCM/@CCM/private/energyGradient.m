function g = energyGradient(c,nx,D,N,T,T_dot,w,W_fcn,dW_fcn) 
    % Compute the gradient of the Riemann Energy under the pseudospectral method
    % taken from Zhao https://github.com/boranzhao/robust_ccm_tube
    persistent c_pre g_pre;
    if isempty(c_pre)  % Adding this may make the gradient calculation inaccurate
        c_pre = zeros(nx,D+1);
        g_pre = zeros(1,(D+1)*nx);
    end
    % gamma = zeros(n,N+1);
    % gamma_s = zeros(n,N+1);
    % for i = 1:n   
    %    gamma(i,:) = c((i-1)*(D+1)+1:i*(D+1),:)'*T;       % gamma(i) is 1*(N+1); the ith elment of gamma on all the (N+1) nodes
    %    gamma_s(i,:) = c((i-1)*(D+1)+1:i*(D+1),:)'*T_dot;
    % end 
    
    c = transpose(reshape(c,D+1,nx)); % The ith row corresponds to the ith element
    if norm(c-c_pre)> 1e-5
        gamma = c*T;
        gamma_s = c*T_dot;
        g = zeros(1,(D+1)*nx);
        dW_dxi = cell(nx,1);
        
        % Vectorized format
        for k = 1:N+1   
            if norm(gamma(3:end,k))> 10
                disp('gamma norm is out of range');
            end
            W_fcn_eval = W_fcn(gamma(:,k));
            M_x_gamma_sk = W_fcn_eval\gamma_s(:,k);
            
            [dW_dxi{:}] = dW_fcn(gamma(:,k));
            for i = 1:nx
                g((i-1)*(D+1)+(1:D+1)) = g((i-1)*(D+1)+(1:D+1))+M_x_gamma_sk'*([zeros(i-1,D+1);T_dot(1:D+1,k)';zeros(nx-i,D+1)]*2-dW_dxi{i}*M_x_gamma_sk*T(1:D+1,k)')*w(k);
            end
        end
        c_pre = c;
        g_pre = g;
    else
        g = g_pre;
    end
end