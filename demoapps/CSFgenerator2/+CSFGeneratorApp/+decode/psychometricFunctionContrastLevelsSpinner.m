function psychometricFunctionContrastLevelsSpinner(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.psychometricFunctionContrastLevelsSpinner.Value = value;
        case 'sliderToValue'
           app.psychometricFunctionParams.contrastLevels = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end