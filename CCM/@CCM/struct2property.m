function obj = struct2property(obj,props,s)
    % Convert input to character cell array
    mustBeText(props)
    props = cellstr(props);
    
    % Split dynamic properties into substrings for get- and setfield
    props = cellfun(@(s) split(s,'.'), props,'UniformOutput',false);
    
    % Loop through struct and assign properties from struct to class
    for i = 1:numel(props)
        try % to set value
            prop = getfield(s, props{i}{:});
            obj = subsasgn(obj, struct('type',repmat({'.'},numel(props{i}),1),'subs',props{i}), prop);
        catch
            id = 'CCM:InvalidFieldorProperty';
            error(id,'%s is not a valid field or property',strjoin(props{i},'.'))
        end
    end
end