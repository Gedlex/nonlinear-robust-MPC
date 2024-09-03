function n_grid = getNgrid(n_vars,grid_pattern,n_tot,options)
    arguments
        n_vars (:,:) {mustBeVector,mustBeInteger}
        grid_pattern (:,:) {mustBeNumeric,mustBePositive}
        n_tot (1,1) {mustBeInteger,mustBePositive} = 1E6
        options.n_grid (:,:) {mustBeVector,mustBeInteger,mustBeNonnegative} = zeros(numel(n_vars),1)
    end
    % Check inputs
    mustBeMember(numel(grid_pattern),numel(n_vars));
    mustBeMember(numel(options.n_grid),numel(n_vars));
    
    % Initialize n_grid
    n_grid = options.n_grid(:);
    
    % Normalize grid pattern
    grid_pattern = grid_pattern./grid_pattern(1);
    
    % Normalize n_tot by already defined grid parameters
    idxs = (n_grid==0);
    if ~all(idxs)
        n_tot = n_tot / prod(n_grid(~idxs).^n_vars(~idxs));
    end
    
    % Compute number of normalized gridding points
    n = (n_tot / prod(grid_pattern(idxs).^n_vars(idxs)))^(1/sum(n_vars(idxs)));
    
    % Compute n_grid
    n_grid(idxs) = n.*grid_pattern(idxs);
    n_grid(n_vars==0) = 1;
    
    % Find all possible combinations
    n_grid_comb = getCombinations([floor(n_grid), ceil(n_grid)]);
    
    % Compute all possible n_tot
    n_tot_comb = prod(n_grid_comb(idxs,:).^n_vars(idxs),1);
    
    % Find minimum
    [~,idx] = min(abs(n_tot_comb - n_tot));
    
    % Return n_grid being closest to n_tot_max
    n_grid = n_grid_comb(:,idx);
end