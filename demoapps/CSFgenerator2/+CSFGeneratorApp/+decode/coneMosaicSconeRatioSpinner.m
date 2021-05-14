function coneMosaicSconeRatioSpinner(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.coneMosaicSconeRatioSpinner.Value = value;
        case 'sliderToValue'
           app.coneMosaicParams.sConeRatio = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end