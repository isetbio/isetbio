function psychometricFunctionSlopeDeltaSpinner(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.psychometricFunctionSlopeDeltaSpinner.Value = value;
        case 'sliderToValue'
           app.psychometricFunctionParams.slopeDelta = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end