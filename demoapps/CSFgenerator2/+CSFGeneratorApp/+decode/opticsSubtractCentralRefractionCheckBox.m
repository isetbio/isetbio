function opticsSubtractCentralRefractionCheckBox(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.opticsSubtractCentralRefractionCheckBox.Value = value;
        case 'sliderToValue'
            app.opticsParams.subtractCentralRefraction = value;
            % Compute new optics
            CSFGeneratorApp.generate.optics(app);
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end