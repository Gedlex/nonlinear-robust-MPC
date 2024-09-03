function vert = poly2vertex(F,b)
    % Check boundedness for constrained states
    vert = [];
    if any((~any(F < 0) | ~any(F > 0)) & any(F))
        warning('Polytope %s is unbounded for constraint state.',inputname(1));
        return
    % Check separability
    elseif any(sum(F ~= 0, 2) > 1)
        warning('Polytope %s cannot be separated into simple 1D boxes.',inputname(1));
        return
    end
    
    % Compute 1D box constraints
    n = size(F,2);
    vert = zeros(n,2);
    for i = 1:n
        % Skip if state is unconstrained
        if ~any(F(:,i))
            continue
        end
        
        % Find nonzero entries
        idx_p = (F(:,i) > 0);
        idx_n = (F(:,i) < 0);
        
        % Find infimum and supremum
        maximum = min(b(idx_p) ./ F(idx_p,i));
        minimum = max(b(idx_n) ./ F(idx_n,i));
        
        % Check that set is nonempty
        if minimum > maximum
            error('Constraint on state %s is empyt set.',num2str(i));
        end
        
        % Add vertices to vert
        vert(i,:) = [minimum, maximum];
    end
end