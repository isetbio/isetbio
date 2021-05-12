function stimulusOrientationSpinner(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.stimulusOrientationSpinner.Value = value;
        case 'sliderToValue'
            app.stimParams.orientationDegs = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end