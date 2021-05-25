function csfSpatialFrequencyMinSpinner(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.csfSpatialFrequencyMinSpinner.Value = value;
        case 'sliderToValue'
           app.csfParams.spatialFrequencyMin = value;
           CSFGeneratorApp.render.csfView(app, 'update');
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
    
    % Update stimulus size at lowest spatial frequency
    stimSizeDegs = CSFGeneratorApp.compute.stimulusSizeAtLowestSpatialFrequency(app);
    app.stimulusSizeAtLowestSpatialFrequencyEditField.Value = sprintf('%2.2f', stimSizeDegs);
end