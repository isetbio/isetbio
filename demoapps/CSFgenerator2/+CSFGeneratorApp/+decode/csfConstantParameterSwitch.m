function csfConstantParameterSwitch(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.csfConstantParameterSwitch.Value = value;
        case 'sliderToValue'
           app.csfParams.constantParameter = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
    
    if (strcmp(value, 'constant size'))
        app.csfNumberOfConstantCyclesSpinner.Enable = 'off';
    else
        app.csfNumberOfConstantCyclesSpinner.Enable = 'on';
    end
    
    % Update stimulus size at lowest spatial frequency
    stimSizeDegs = CSFGeneratorApp.compute.stimulusSizeAtLowestSpatialFrequency(app);
    app.stimulusSizeAtLowestSpatialFrequencyEditField.Value = sprintf('%2.2f', stimSizeDegs);
end