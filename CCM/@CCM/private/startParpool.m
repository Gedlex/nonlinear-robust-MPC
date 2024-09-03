function poolobj = startParpool(profile)
    arguments
        profile (1,1) string = parallel.defaultClusterProfile
    end
    poolobj = gcp("nocreate"); % If no pool, do not create new one.
    if isempty(poolobj)
        poolobj = parpool(profile);
    elseif ~strcmp(poolobj.Cluster.Profile,profile)
        delete(poolobj)
        poolobj = parpool(profile);
    end
end