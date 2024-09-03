function [f_sym, sym_vars] = sdpvar2syms(varargin)
    % Copy sdpvar to function workspace
    for i = 1:nargin
        if ~isa(varargin{i},'sdpvar')
            error('All input arguments must be sdpvar objects')
        elseif i==1
            continue
        end
        var_name = inputname(i);
        eval([var_name,'=varargin{i};'])
    end
    
    % Convert sdpvar to string
    svar = sdisplay(varargin{1});
    
    if any(contains(svar,'internal'))
        error(['Not enough input arguments. Provide all sdpvar used in ',inputname(1),' as an input.'])
    end
    
    % Create symbolic variables from input
    clearvars('-except','varargin','svar')
    sym_vars = cell(nargin-1,1);
    for i = 2:nargin
        var_name = inputname(i);
        eval([var_name,'=sym(var_name, size(varargin{i}));'])
        sym_vars{i-1} = eval(var_name);
    end
    
    % Convert string to symbolic expression
    f_sym = sym('f_sim', size(svar));
    for i = 1:numel(f_sym)
        f_sym(i) = eval(svar{i});
    end
end