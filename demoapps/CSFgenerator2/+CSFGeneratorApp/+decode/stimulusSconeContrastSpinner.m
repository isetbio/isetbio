function stimulusSconeContrastSpinner(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.stimulusSconeContrastSpinner.Value = value;
        case 'sliderToValue'
            app.stimParams.SconeContrast = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end