function psychometricFunctionEstimationMethodSwitch(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.psychometricFunctionEstimationMethodSwitch.Value = value;
        case 'sliderToValue'
           app.psychometricFunctionParams.estimationMethod = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
    
    if (contains(value, 'constant stimuli')) 
        app.psychometricFunctionSlopeMinSpinner.Enable = 'off';
        app.psychometricFunctionSlopeMaxSpinner.Enable = 'off';
        app.psychometricFunctionSlopeDeltaSpinner.Enable = 'off';
    else
        app.psychometricFunctionSlopeMinSpinner.Enable = 'on';
        app.psychometricFunctionSlopeMaxSpinner.Enable = 'on';
        app.psychometricFunctionSlopeDeltaSpinner.Enable = 'on';
    end
end