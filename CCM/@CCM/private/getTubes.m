function tubes = getTubes(F,tightenings)
    tubes = [];
    if any(~any(F < 0) | ~any(F > 0)) || any(sum(F ~= 0, 2) > 1)
        warning('Can not compute simple 1D tubes from tightenings.');
        return
    end
    
    % Compute tubes
    n = size(F,2);
    m = size(tightenings,2);
    tubes = zeros(n,m);
    for i = 1:n        
        % Find nonzero entries
        idx = (F(:,i) ~= 0);
        
        % Compute maximum tightening
        tubes(i,:) = max(abs(tightenings(idx,:) ./ F(idx,i)),[],1);
    end
end