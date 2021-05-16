function psychometricFunctionClassifierTypeDropDown(app, direction, value)
    
    switch direction
        case 'valueToSlider'
            app.psychometricFunctionClassifierTypeDropDown.Value = value;
        case 'sliderToValue'
            app.psychometricFunctionParams.classifierType = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
    
    if (contains(value, 'ideal observer'))
        app.psychometricFunctionClassifierTrainingTrialsSpinner.Enable = 'off';
    else
        app.psychometricFunctionClassifierTrainingTrialsSpinner.Enable = 'on';
    end
end

