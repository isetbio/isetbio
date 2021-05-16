function psychometricFunctionClassifierTestTrialsSpinner(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.psychometricFunctionClassifierTestTrialsSpinner.Value = value;
        case 'sliderToValue'
           app.psychometricFunctionParams.testTrials = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end