function psychometricFunctionLog10ContrastMaxSpinner(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.psychometricFunctionLog10ContrastMaxSpinner.Value = value;
        case 'sliderToValue'
           app.psychometricFunctionParams.log10ContrastMax = value;
           CSFGeneratorApp.render.psychometricFunctionView(app, 'update');
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end