function allWarnings(state)

    return;
    
    if (nargin == 0) || ((nargin == 1) && (~ischar(state)))
        printWarningStates();
        return;
    end
    switch (state)
        case 'default'
            setDefaultState();
        case 'on'
            setAll('on');
        case 'off'
            setAll('off');
        otherwise
            fprintf('Unrecognized state value (''%s''). Choose from: {''default'', ''on'', ''off''}.', state);
    end
    printWarningStates();
end

function printWarningStates()
    warningLabels = coneMosaic.warnings.keys;
    for k = 1:numel(warningLabels)
        warning('query', warningLabels{k});
    end
end

function setAll(state)
    warningLabels = coneMosaic.warnings.keys;
    for k = 1:numel(warningLabels)
        warning(state, warningLabels{k});
    end
end

function setDefaultState()
    warningLabels = coneMosaic.warnings.keys;
    warningStates = coneMosaic.warnings.values;
    for k = 1:numel(warningLabels)
        warning(warningStates{k}, warningLabels{k});
    end
end

