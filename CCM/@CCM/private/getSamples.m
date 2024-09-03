function samples = getSamples(F,b,n_grid,keep)
    arguments
        F (:,:) {mustBeNumeric}
        b (:,:) {mustBeNumeric,mustBeVector}
        n_grid (1,1) {mustBeInteger,mustBePositive}
        keep (:,:) logical {mustBeVector} = true(size(F,2),1)
    end
    
    % Convert polytope to vertices
    vertices = poly2vertex(F,b);
    n = size(F,2);
    
    % Check if conversion succeeded
    if ~isempty(vertices)
        % Initialize zero vector
        samples = repmat({0},n,1);
        
        % Grid over vertices in keep
        for k = 1:n
            if keep(k)
                samples{k} = linspace(vertices(k,1),vertices(k,2),n_grid);
            end
        end
        
        % Permute vector to get all grid combinations
        samples = getCombinations(samples{:});
    else
        % Inform user of possible inaccuracies
        warning(['Polytope can`t be separated and must be sampled at random. ' ...
                 'This might be less accurate than gridding. ' ...
                 'Check the result carefully or consider defining polytopes that can be gridded.'])

        % Fix random seed for reproducable results
        rng(0,'twister');
        
        % Get unconstrained states
        constr = any(F,1)';
        
        % Compute number of samples (to be equivalent to 2*gridding)
        n_sample = 2*n_grid^(sum(keep(:) & constr));
        
        % Initialize zero vector
        samples = zeros(n,n_sample);
        
        % Sample from constrained polyhedron only
        samples(constr,:) = cprnd(n_sample, F(:,constr), b)';
    end
end