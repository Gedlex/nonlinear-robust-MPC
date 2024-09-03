function s = property2struct(obj,props,s)
    % Convert input to character cell array
    mustBeText(props)
    props = cellstr(props);
    
    % Split dynamic properties into substrings for get- and setfield
    props = cellfun(@(s) split(s,'.'), props,'UniformOutput',false);
    
    % Allocate output struct
    if ~exist('s','var')
        s = struct();
    end
    
    % Loop through properties and append them to the struct
    for i = 1:numel(props)
        try % to set value
            prop = getfield(obj, props{i}{:});
            s = setfield(s, props{i}{:}, prop);
        catch
            id = 'CCM:InvalidFieldorProperty';
            error(id,'%s is not a valid field or property',strjoin(props{i},'.'))
        end
    end
end