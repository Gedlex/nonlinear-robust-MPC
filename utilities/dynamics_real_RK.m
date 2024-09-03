function [x_new, u1] = dynamics_real_RK(f,x,kappa,z,v,d,h)
    arguments
        f (1,1) {mustBeA(f,'function_handle')}
        x (:,:) {mustBeVector}
        kappa (1,1) {mustBeA(kappa,'function_handle')}
        z (:,:) {mustBeNumeric}
        v (:,:) {mustBeNumeric}
        d (:,:) {mustBeNumeric,mustBeVector}
        h (1,1) {mustBePositive} = 0.1
    end
    % Simulates the true dynamics using ruku4
    % f: function handle to dynamics
    % x: true state at time t_k
    % kappa: function handle to feedback control law
    % z: nominal states stacked [z_{t_k}, z_{t_k + h/2}, z_{t_k + h}]  
    % v: nominal states stacked [v_{t_k}, v_{t_k + h/2}, v_{t_k + h}]  
    % d: disturbance realizations stacked [d_{t_k}, d_{t_k + h/2}, d_{t_k + h}] 
    % h: ruku4 step size
    % returns x_new: true state at t_{k+1}
    % returns u1: applied input at time t_k
    
    u1 = kappa(x,z(:,1),v(:,1));
    k1=f(x,u1,d(:,1));
    u2 = kappa(x+h/2*k1,z(:,2),v(:,2));
    k2=f(x+h/2*k1,u2,d(:,2));
    u3 = kappa(x+h/2*k2,z(:,2),v(:,2));
    k3=f(x+h/2*k2,u3,d(:,2));
    u4 = kappa(x+h*k3,z(:,3),v(:,3));
    k4=f(x+h*k3,u4,d(:,3));
    x_new=x+h/6*(k1+2*k2+2*k3+k4);
end