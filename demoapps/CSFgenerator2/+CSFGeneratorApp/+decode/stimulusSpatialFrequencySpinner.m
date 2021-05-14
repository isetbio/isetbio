function stimulusSpatialFrequencySpinner(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.stimulusSpatialFrequencySpinner.Value = value;
        case 'sliderToValue'
            app.stimParams.spatialFrequencyCPD = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end