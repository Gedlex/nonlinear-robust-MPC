function field_name = findMatchingField(s1,s2,identifier,exactMatch)
    % Get struct fieldnames and allocate return value
	fields = fieldnames(s1);
    field_name = '';

    % Boolean for non-existing fields
    if ~exist('exactMatch','var')
        exactMatch = false;
    end
    
    % Loop through struct and search for matching entry
    for i = 1:numel(fields)
        % Check only field names that start with the correct identifier
        if ~startsWith(fields{i},identifier)
            continue
        end
		% Check if structs match
        if ~structcmp(s2,s1.(fields{i}),exactMatch)
			continue
        end
		% Return field name found
		field_name = fields{i};
		break
    end
end