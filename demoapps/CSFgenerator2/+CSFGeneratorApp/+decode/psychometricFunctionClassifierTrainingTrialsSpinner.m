function psychometricFunctionClassifierTrainingTrialsSpinner(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.psychometricFunctionClassifierTrainingTrialsSpinner.Value = value;
        case 'sliderToValue'
           app.psychometricFunctionParams.trainingTrials = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end