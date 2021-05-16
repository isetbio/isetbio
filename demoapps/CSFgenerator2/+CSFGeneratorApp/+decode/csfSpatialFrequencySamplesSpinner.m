function csfSpatialFrequencySamplesSpinner(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.csfSpatialFrequencySamplesSpinner.Value = value;
        case 'sliderToValue'
           app.csfParams.spatialFrequencySamples = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end