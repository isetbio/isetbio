function psychometricFunctionSlopeMaxSpinner(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.psychometricFunctionSlopeMaxSpinner.Value = value;
        case 'sliderToValue'
           app.psychometricFunctionParams.slopeMax = value;
           CSFGeneratorApp.render.psychometricFunctionView(app, 'update');
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end