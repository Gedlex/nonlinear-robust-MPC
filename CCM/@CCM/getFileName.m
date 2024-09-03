function file_name = getFileName(obj)
    % Extract class directory
    [path, ~, ~] = fileparts(mfilename('fullpath'));
    
    % Return file name
    file_name = fullfile(path,[class(obj.sys),'.mat']);
end