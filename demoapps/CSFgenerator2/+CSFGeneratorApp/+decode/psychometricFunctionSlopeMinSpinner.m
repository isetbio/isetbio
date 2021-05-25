function psychometricFunctionSlopeMinSpinner(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.psychometricFunctionSlopeMinSpinner.Value = value;
        case 'sliderToValue'
           app.psychometricFunctionParams.slopeMin = value;
           CSFGeneratorApp.render.psychometricFunctionView(app, 'update');
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end