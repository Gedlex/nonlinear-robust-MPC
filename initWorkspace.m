function initWorkspace()
    % Extract directory
    [path, ~, ~] = fileparts(mfilename('fullpath'));
    
    % Add directory and all subfolders to path
    addpath(genpath(path));
    cd(path)
end