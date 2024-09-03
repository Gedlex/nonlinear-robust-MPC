function saveVars(obj,name,vars,props)    
	% Check if file already exists
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
        props_struct = obj.property2struct(props);
        
        % Find field with prefix name
        field_name = findMatchingField(data,props_struct,name);
        
        % Check if field was found
        props = [vars, props];
        if isempty(field_name)
            % Create unique field name
            field_name = getUniqueFieldName(data,name);
            
            % Create new struct from properties
            s.(field_name) = obj.property2struct(props);
        else
            % Add properties to existing struct
            s.(field_name) = obj.property2struct(props,data.(field_name)); 
        end
        
        % Wait for lock
        while isfile(lock)
            pause(rand)
        end
        
        % Create lock
        fid = fopen(lock,'w');
        fclose(fid);
        
        % Save offline computations
        save(filepath,'-struct','s',field_name,'-append');
        
        % Delete lock
        delete(lock)
    else
        % Create lock
        fid = fopen(lock,'w');
        fclose(fid);
        
        % Save offline computations
        field_name = [name,'_0'];
        props = [vars, props];
        s.(field_name) = obj.property2struct(props);
        save(filepath,'-struct','s',field_name);
        
        % Delete lock
        delete(lock)
    end
end