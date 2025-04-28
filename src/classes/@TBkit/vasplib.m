function vasplibobj =vasplib(propArgs)
    arguments
        propArgs.?vasplib;
    end
    Fieldnames = fieldnames(propArgs);
    if isempty(Fieldnames)
    else
    for iFieldname = 1:numel(Fieldnames)
        vasplibobj.(Fieldnames{iFieldname}) = propArgs.(Fieldnames{iFieldname});
    end
    end
    vasplibobj = vasplibobj.vasplib_init();
end
