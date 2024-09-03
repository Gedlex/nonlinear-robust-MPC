function Theta = compute_non_falsified(meas_x,meas_u,meas_dx,d_lim,eps_lim,m_nom)
    % Compute non-falsified set for set-membership estimation
    % Input: meas_x (measured state), meas_u (measured input), meas_dx (noisy measured state derivative), d_lim (bound maximal disturbance), eps_lim (bound on measurement noise), params (model parameters)
    % Output: Theta (non-falsified set, written in terms of lower and upper bound)
    % Code is tailored for specific scalar case
    tmp1 = [];
    for eps = linspace(-eps_lim,eps_lim,2)        
        for d = linspace(-d_lim,d_lim,2)
            tmp1 = [tmp1; (meas_dx(5,:) + meas_x(4,:).*meas_x(6,:) + 9.81*cos(meas_x(3,:)) ...
                + sin(meas_x(3,:))*d - eps)./([1 1]*meas_u)];
        end
    end
    
    tmp1_min = min(tmp1,[],1);
    tmp1_max = max(tmp1,[],1);
    
    Theta = [max(tmp1_min); min(tmp1_max)] - 1/m_nom;
end