function progbar(progress,options)
    arguments
        progress (1,1) {mustBeInRange(progress,0,100)} = 0
        options.increment (1,1) {mustBeNonnegative} = 0
        options.barLength (1,1) {mustBeInteger} = 0
        options.reset (1,1) logical = 0
        options.prefix (1,1) string = ""
    end
    persistent storedProgress prefix barLength
    
    % Initialize persistent variables
    if isempty(storedProgress) || options.reset
        storedProgress = 0;
        prefix = 'Progress: ';
        barLength = 25;
        if options.reset; return; end
    end
    
    % Update by increment
    if options.increment && ~progress
        progress = min(storedProgress + options.increment, 100);
    end
    
    % Generate reverse string
    strLength = 0;
    if 0 < progress && storedProgress <= progress
        strLength = strlength(prefix) + barLength + 8;
    end
    reverseStr = repmat(sprintf('\b'), 1, strLength);
    
    % Change persistent variables
    if options.prefix ~= ""
        prefix = options.prefix;
    end
    if options.barLength
        barLength = options.barLength;
    end
    
    % Build bar
    n_bars = round(progress*barLength/100);
    bar = sprintf('[%s%s]',repmat('#',1,n_bars),repmat('-',1,barLength-n_bars));
    
    % Print progress
    msg = sprintf('%s%s %3u%% ',prefix,bar,round(progress));
    fprintf([reverseStr,'%s'],msg);

    % Check for end
    if (100-progress < 1E-5) && (100-storedProgress > 1E-5)
        fprintf('\n');
    end

    % Update stored progress
    storedProgress = progress;
end