function roiFieldOfViewSlider(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.roiFieldOfViewSlider.Value = value;
        case 'sliderToValue'
            app.roiParams.fieldOfViewDegs = value;
            if (value == 0)
                app.roiParams.fieldOfViewDegs = 0.1;
            end
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end

    % Update the text
    if (app.roiParams.fieldOfViewDegs < 1.0)
        app.fieldofviewdegsLabel.Text = sprintf('field of view (degs): %2.2f x %2.2f', app.roiParams.fieldOfViewDegs, 0.5*app.roiParams.fieldOfViewDegs);
    elseif (app.roiParams.fieldOfViewDegs < 10)
        app.fieldofviewdegsLabel.Text = sprintf('field of view (degs): %2.1f x %2.1f', app.roiParams.fieldOfViewDegs, 0.5*app.roiParams.fieldOfViewDegs);
    else
        app.fieldofviewdegsLabel.Text = sprintf('field of view (degs): %2.0f x %2.0f', app.roiParams.fieldOfViewDegs, 0.5*app.roiParams.fieldOfViewDegs);
    end

    % Set the associated coneMosaicParam
    app.coneMosaicParams.sizeDegs = app.roiParams.fieldOfViewDegs * [1 0.5];
end
