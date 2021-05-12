function opticsPupilSizeSpinner(app, direction, value)
    
    switch direction
        case 'valueToSlider'
            app.opticsPupilSizeSpinner.Value = value;
        case 'sliderToValue'
            app.opticsParams.pupilDiameterMM = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end