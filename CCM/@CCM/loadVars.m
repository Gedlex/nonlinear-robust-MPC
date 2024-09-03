function [obj, success] = loadVars(obj,name,vars,props)
    % Initialize return
    success = false;
    
    % Check if file exists
    filepath = obj.getFileName();
    lock = [filepath,'.lock'];
    if isfile(filepath)       
        % Wait for lock
        while isfile(lock)
            pause(rand)
        end
        
        % Create lock
        fid = fopen(lock,'w');
        fclose(fid);
        
        % Load structured array
        data = load(filepath);
        
        % Delete lock
        delete(lock)
        
        % Convert properties to structured array
        props = obj.property2struct(props);
        
        % Find field with prefix name that matches all properties in props exactly
        field_name = findMatchingField(data,props,name,true);
        
        % Load vars from file
        if ~isempty(field_name)
            try % to load variables from file
                obj = obj.struct2property(vars,data.(field_name));
                success = true;
            catch
                return
            end
        end
    end
end