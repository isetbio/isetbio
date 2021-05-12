function stimulusResolutionSpinner(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.stimulusResolutionSpinner.Value = value;
        case 'sliderToValue'
            app.stimParams.resolutionPixels = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end