function opticsWavefrontSpatialSamplesSpinner(app, direction, value)
    
    switch direction
        case 'valueToSlider'
            app.opticsWavefrontSpatialSamplesSpinner.Value = value;
        case 'sliderToValue'
            app.opticsParams.wavefrontSpatialSamples = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end