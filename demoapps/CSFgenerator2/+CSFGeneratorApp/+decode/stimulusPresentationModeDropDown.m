function stimulusPresentationModeDropDown(app, direction, value)
    
    switch direction
        case 'valueToSlider'
            app.stimulusPresentationModeDropDown.Value = value;
        case 'sliderToValue'
            app.stimParams.presentationMode = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end

