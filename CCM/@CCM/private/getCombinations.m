function comb = getCombinations(varargin)
    % Convert stacked vectors into cell arrays
    cellvec = cellfun(@(X) mat2cell(X, ones(1,size(X,1))), varargin, 'UniformOutput', false)';
    cellvec = vertcat(cellvec{:}); % flatten nested cells
    
    % Remove duplicates from vector
    cellvec = cellfun(@(X) unique(X), cellvec, 'UniformOutput', false);
    
    % Compute n-dimensional grid
    out = cell(size(cellvec,1),1);
    [out{:}] = ndgrid(cellvec{:});
    
    % Convert ndgrid into permuation matrix
    out = cellfun(@(X) X(:)', out, 'UniformOutput', false);
    comb = cell2mat(out);
end