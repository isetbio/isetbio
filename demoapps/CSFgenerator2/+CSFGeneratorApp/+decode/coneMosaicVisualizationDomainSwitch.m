function coneMosaicVisualizationDomainSwitch(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.coneMosaicVisualizationDomainSwitch.Value = value;
        case 'sliderToValue'
            app.viewModes.coneMosaicVisualizationDomain = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end
