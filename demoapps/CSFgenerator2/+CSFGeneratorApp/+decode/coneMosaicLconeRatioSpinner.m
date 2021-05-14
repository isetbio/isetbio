function coneMosaicLconeRatioSpinner(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.coneMosaicLconeRatioSpinner.Value = value;
        case 'sliderToValue'
           app.coneMosaicParams.lConeRatio = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end