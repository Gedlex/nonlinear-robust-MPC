function E = RiemannEnergy(c,nx,D,N,T,T_dot,w,W_fcn)
    % Compute the Riemann Energy using the pseudospectral method
    % taken from Zhao https://github.com/boranzhao/robust_ccm_tube
    
    % Vectorized format to improve computational efficiency
    c = transpose(reshape(c,D+1,nx)); % the ith row corresponds to the ith element
    
    gamma = c*T;
    gamma_s = c*T_dot;
    E = 0;
    for k=1:N+1  
        tmp = gamma_s(:,k)'*(W_fcn(gamma(:,k))\gamma_s(:,k))*w(k);
        E = E + tmp;
    end
end