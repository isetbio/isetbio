function coneMosaicTritanopicRadiusSpinner(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.coneMosaicTritanopicRadiusSpinner.Value = value;
        case 'sliderToValue'
           app.coneMosaicParams.tritanopicRadiusDegs = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end