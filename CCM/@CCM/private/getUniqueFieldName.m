function unique_field_name = getUniqueFieldName(s,field_name)
    suffix = 0;
    % Test suffix until unique
    unique_field_name = sprintf('%s_%d',field_name,suffix);
    while isfield(s, unique_field_name)
        unique_field_name = sprintf('%s_%d',field_name,suffix);
        suffix = suffix + 1;
    end
end