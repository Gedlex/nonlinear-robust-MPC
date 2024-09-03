function TF = structcmp(s1,s2,exactMatch)
    % Get struct fieldnames and allocate return value
    fields = fieldnames(s1);
    TF = false;
    
    % Compare structs against each other
    for i = 1:numel(fields)
        % Check if field exists
        if ~isfield(s2, fields{i})
            if exactMatch % return false
                return
            else % compare against rest
                continue
            end
        end
        % Recursively call structcmp to compare nested structs
        if isstruct(s1.(fields{i})) && isstruct(s2.(fields{i}))
            if ~structcmp(s1.(fields{i}),s2.(fields{i}),exactMatch)
                return
            end
        % Compare other data (does not work for anonymouse functions)
        elseif ~isequal(s1.(fields{i}), s2.(fields{i}))
            return
        end
    end
    TF = true;
end