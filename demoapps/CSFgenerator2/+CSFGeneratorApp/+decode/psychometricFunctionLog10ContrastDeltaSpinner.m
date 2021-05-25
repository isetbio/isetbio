function psychometricFunctionLog10ContrastDeltaSpinner(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.psychometricFunctionLog10ContrastDeltaSpinner.Value = value;
        case 'sliderToValue'
           app.psychometricFunctionParams.log10ContrastDelta = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end