function psychometricFunctionLog10ContrastMinSpinner(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.psychometricFunctionLog10ContrastMinSpinner.Value = value;
        case 'sliderToValue'
           app.psychometricFunctionParams.log10ContrastMin = value;
           CSFGeneratorApp.render.psychometricFunctionView(app, 'update');
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end