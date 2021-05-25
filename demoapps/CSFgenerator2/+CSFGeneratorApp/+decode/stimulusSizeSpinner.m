function stimulusSizeSpinner(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.stimulusSizeSpinner.Value = value;
        case 'sliderToValue'
            app.stimParams.sizeDegs = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
    
    % Update stimulus size at lowest spatial frequency
    stimSizeDegs = CSFGeneratorApp.compute.stimulusSizeAtLowestSpatialFrequency(app);
    app.stimulusSizeAtLowestSpatialFrequencyEditField.Value = sprintf('%2.2f', stimSizeDegs);
end