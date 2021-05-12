function stimulusSizeSpinner(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.stimulusSizeSpinner.Value = value;
        case 'sliderToValue'
            app.stimParams.sizeDegs = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end