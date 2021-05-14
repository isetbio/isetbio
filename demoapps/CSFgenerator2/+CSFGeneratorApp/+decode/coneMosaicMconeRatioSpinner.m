function coneMosaicMconeRatioSpinner(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.coneMosaicMconeRatioSpinner.Value = value;
        case 'sliderToValue'
           app.coneMosaicParams.mConeRatio = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end