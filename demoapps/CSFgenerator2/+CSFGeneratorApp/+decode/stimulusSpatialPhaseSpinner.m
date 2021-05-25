function stimulusSpatialPhaseSpinner(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.stimulusSpatialPhaseSpinner.Value = value;
        case 'sliderToValue'
            app.stimParams.spatialPhaseDegs = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end