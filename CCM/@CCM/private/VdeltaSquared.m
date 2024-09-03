function Vd = VdeltaSquared(N,gamma,gamma_s,w,W_fcn)
    % Integrate the Riemannian Energy E = V_delta^2 (see Eqn. (8) and App. A)
    % based on code from Zhao https://github.com/boranzhao/robust_ccm_tube
    Vd = 0;
    for k=1:N+1  
        tmp = gamma_s(:,k)'*(W_fcn(gamma(:,k))\gamma_s(:,k))*w(k);
        Vd = Vd + tmp ;
    end
end