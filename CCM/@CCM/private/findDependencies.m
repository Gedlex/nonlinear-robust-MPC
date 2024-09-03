function TF = findDependencies(fcn,varargin)
    import casadi.*
    % Check function
    n = numel(varargin);
    if isa(fcn,'function_handle') || isa(fcn,'casadi.Function')
        % Check input
        cellfun(@(X) mustBeNumeric(X), varargin);
        
        % Create symbolic variables
        sym_vars = cell(n,1);
        for i = 1:n
            % Get size from input
            inp = varargin{i};
            if isscalar(inp)
                dim = [inp,1];
            elseif numel(inp) == 2
                dim = inp(:);
            else
                dim = size(inp);
            end
            
            % Create symbolic variable
            name = sprintf('x%u',i);
            sym_vars{i} = SX.sym(name,dim);
        end
        
        % Transform function to symbolic expression
        sym_fcn = fcn(sym_vars{:});
    else
        % Check input
        mustBeA(fcn,["casadi.SX","casadi.MX","casadi.DM"]);
        cellfun(@(X) mustBeA(X,["casadi.SX","casadi.MX","casadi.DM"]), varargin);
        
        % Pass inputs on
        sym_fcn = fcn;
        sym_vars = varargin;
    end
    
    % Determine symbolic variables used in sym_fcn
    TF = cell(n,1);
    for i = 1:n
        TF{i} = false(size(sym_vars{i}));
        for j = 1:numel(sym_vars{i})
            TF{i}(j) = depends_on(sym_fcn, sym_vars{i}(j));
        end
    end
end