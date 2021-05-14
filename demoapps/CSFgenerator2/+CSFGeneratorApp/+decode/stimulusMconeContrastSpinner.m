function stimulusMconeContrastSpinner(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.stimulusMconeContrastSpinner.Value = value;
        case 'sliderToValue'
            app.stimParams.MconeContrast = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end