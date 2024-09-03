function box = poly2box(F_x,b_x,F_u,b_u,F_d,b_d,theta_v,options)
    arguments
        F_x (:,:) {mustBeNumeric}
        b_x (:,:) {mustBeNumeric,mustBeVector}
        F_u (:,:) {mustBeNumeric}
        b_u (:,:) {mustBeNumeric,mustBeVector}
        F_d (:,:) {mustBeNumeric}
        b_d (:,:) {mustBeNumeric,mustBeVector}
        theta_v (:,2) {mustBeNumeric}
        options.normalize (1,1) logical = true
        options.scaling (:,:) {mustBeNumeric,mustBeVector} = ones(size(F_x,2)+size(F_u,2)+size(F_d,2)+size(theta_v,1),1)
    end
    % Check polytopes fo boundedness and separability
    polytopes = {F_x,F_u,F_d};
    for i = 1:numel(polytopes)
        currentPolytope = polytopes{i};

        % Check boundedness for constrained states
        if any((~any(currentPolytope < 0) | ~any(currentPolytope > 0)) & any(currentPolytope))
            error('Polytope %s is unbounded for constraint state.',inputname(2*i-1));

        % Check separability
        elseif any(sum(currentPolytope ~= 0, 2) > 1)
            error('Polytope %s cannot be separated into simple 1D ellipsoidal boxes.',inputname(2*i-1));
        end
    end

    % Get dimensions
    nx = size(F_x,2);
    nu = size(F_u,2);
    nw = size(F_d,2);
    np = size(theta_v,1);
    
    % Create symbolic variables
    x_sym = sym('x',[nx,1]);
    u_sym = sym('u',[nu,1]);
    d_sym = sym('d',[nw,1]);
    theta_sym = sym('theta',[np,1]);
    sym_vars = [x_sym;u_sym;theta_sym;d_sym];
    
    % Concatenate polyhedrons
    F_tot = blkdiag(F_x,F_u,[eye(np);          -eye(np)], F_d);
    b_tot = [b_x(:); b_u(:);[theta_v(:,2);-theta_v(:,1)]; b_d];
    
    % Compute box constraint
    box=[];
    for i = 1:size(F_tot,2)
        % Skip if state is unconstrained (blank column)
        if ~any(F_tot(:,i))
            continue
        end
        
        % Find nonzero entries
        idx_p = (F_tot(:,i) > 0);
        idx_n = (F_tot(:,i) < 0);
        
        % Find infimum and supremum
        maximum = min(b_tot(idx_p) ./ F_tot(idx_p,i));
        minimum = max(b_tot(idx_n) ./ F_tot(idx_n,i));
        
        % Check that set is nonempty
        if minimum > maximum
            error('Constraint on %s is empyt set.',string(sym_vars(i)));
        end
        
        % Compute mean and ellipsoidal bound
        shift = mean([minimum, maximum]);
        bound = maximum - shift;
        
        % Compute normalizer
        normalizer = 1;
        if options.normalize && 1E-10 < bound^2
            normalizer = bound^2;
        end
        
        % Add constraint to box contraint
        box = [box; options.scaling(i)*(bound^2 - (sym_vars(i) - shift)^2) / normalizer];
    end
    
    % Convert to variable-precision arithmetic representation
    box = vpa(box);
    
    % Convert to matlab function
    box = matlabFunction(box,'Vars',{x_sym,u_sym,theta_sym,d_sym});
end