function fieldOfViewSlider(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.visualFieldFieldOfViewSlider.Value = value;
        case 'sliderToValue'
            app.visualFieldParams.fieldOfViewDegs = value;
            if (value == 0)
                app.visualFieldParams.fieldOfViewDegs = 0.1;
            end
        otherwise
            error('Unknown FieldOfViewSlider.direction; ''%s''.'\n', direction);
    end

    % Update the text
    if (app.visualFieldParams.fieldOfViewDegs < 1.0)
        app.fieldofviewdegsLabel.Text = sprintf('field of view (degs): %2.2f', app.visualFieldParams.fieldOfViewDegs);
    elseif (app.visualFieldParams.fieldOfViewDegs < 10)
        app.fieldofviewdegsLabel.Text = sprintf('field of view (degs): %2.1f', app.visualFieldParams.fieldOfViewDegs);
    else
        app.fieldofviewdegsLabel.Text = sprintf('field of view (degs): %2.0f', app.visualFieldParams.fieldOfViewDegs);
    end

    % Set the associated coneMosaicParam
    app.coneMosaicParams.sizeDegs = app.visualFieldParams.fieldOfViewDegs * [1 1];
end
