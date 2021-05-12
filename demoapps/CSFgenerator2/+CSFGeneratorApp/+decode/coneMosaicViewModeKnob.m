function coneMosaicViewModeKnob(app, direction, value)
    
    switch direction
        case 'valueToSlider'
            app.coneMosaicViewModeKnob.Value = value;
        case 'sliderToValue'
            app.viewModes.coneMosaic = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
    
end
