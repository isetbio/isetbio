function stimulusLconeContrastSpinner(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.stimulusLconeContrastSpinner.Value = value;
        case 'sliderToValue'
            app.stimParams.LconeContrast = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end