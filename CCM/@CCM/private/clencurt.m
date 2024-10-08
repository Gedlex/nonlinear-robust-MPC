function [x,w] = clencurt(N)
    % CLENCURT  nodes x (Chebyshev points) and weights w for Clenshaw-Curtis
    %           quadrature (taken from Zhao https://github.com/boranzhao/robust_ccm_tube)
    
    theta = pi*(0:N)'/N; x = cos(theta);
    
    w = zeros(1,N+1); ii = 2:N; v = ones(N-1,1);
    if mod(N,2)==0 
    w(1) = 1/(N^2-1); w(N+1) = w(1);
    for k=1:N/2-1, v = v - 2*cos(2*k*theta(ii))/(4*k^2-1); end
    v = v - cos(N*theta(ii))/(N^2-1);
    else
    w(1) = 1/N^2; w(N+1) = w(1);
    for k=1:(N-1)/2, v = v - 2*cos(2*k*theta(ii))/(4*k^2-1); end
    end
    w(ii) = 2*v/N;
    
    x = flip(x); w = flip(w); x=x';
    
    x=(x+1)/2;  % Convert the integration range from [-1, 1] to [0,1]
    w = w/2;    % Accommodate the above change
end