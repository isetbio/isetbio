function csfSpatialFrequencyMaxSpinner(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.csfSpatialFrequencyMaxSpinner.Value = value;
        case 'sliderToValue'
           app.csfParams.spatialFrequencyMax = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end