function u = kappa(obj,x,z,v)
    % Compute geodesic 
    [gamma,gamma_s] = computeGeodesic(x,z,obj.geod);
    
    % Integrate feedback control law kappa
    u = v;
    for k=1:obj.geod.N+1 
        u = u + obj.geod.w_cheby(k)*(obj.Y_fcn(gamma(:,k))*(obj.W_fcn(gamma(:,k))\gamma_s(:,k)));
    end
end