function ret = w_tilde(obj,x,u,theta_bar,theta,d,W_chol_x)
    % Compute || G(u)*(theta-theta_bar) + E(x)*d ||_{M(z)}.
    % v'*M(z)*v = v'*M^1/2(z)'*M^1/2(z)*v, where M^1/2 can be any 
    % matrix square root. We use M^1/2(z) = chol(W(z)).
    vec = obj.sys.G_fcn(x,u)*(theta - theta_bar) + obj.sys.E_fcn(x)*d;
    ret = norm(vec'/W_chol_x);
end