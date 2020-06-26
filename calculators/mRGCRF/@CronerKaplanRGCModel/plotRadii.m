function plotRadii(theAxes, d, model, pointSize, color, displayYLabel, theLabel)

    d = d('size');
    
    % Face color
    c = color + 0.5*[1 1 1];
    c(c>1)=1;
    
    hold(theAxes, 'on');
    
    % Plot the digitized data
    if (strcmp(theLabel, 'center')) || (strcmp(theLabel, 'surround'))
        scatter(theAxes, d.eccDegs, d.radiusDegs, pointSize, 'MarkerFaceColor', c, 'MarkerEdgeColor', color);
    elseif (strcmp(theLabel, 'retinal center'))
        scatter(theAxes, d.retinalEccDegs, d.retinalRadiusDegs, pointSize, 'MarkerFaceColor', c, 'MarkerEdgeColor', color);
    else
        error('Label must b either ''center'', ''surround'', or ''retinal center''.');
    end
    
    if (isfield(d, 'eccDegsTable')) &&  (~(strcmp(theLabel, 'retinal center')))
        scatter(theAxes, d.eccDegsTable, d.radiusDegsMedianTable, 's', 'MarkerFaceColor', [0.2 0.2 0.2], 'MarkerEdgeColor', [0.2 0.2 0.2]);
    end
    
    % Plot the fitted model
    plot(theAxes, model.eccDegs, model.function(model.params, model.eccDegs), 'k-');
    
    % Add text with model equation
    t=['$\displaystyle \rho = ', sprintf('%2.3f',model.params(1)),' \times {\epsilon}^{', sprintf('%2.3f',model.params(2)), '}$'];
    
    theTextHandle = text(theAxes, 0.2, 7 , t, 'Interpreter', 'latex');
    set(theTextHandle,'FontSize', 20, 'Color', [0.3 0.3 0.3], 'BackgroundColor', [1 1 1]);

    % Finish plot
    set(theAxes, 'XLim', [0.1 100], 'YLim', [0.01 10]);
    set(theAxes, 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30 100], 'YTick', [0.01 0.03 0.1 0.3 1 3 10]);
    set(theAxes, 'XTickLabel', {'0.01' '0.03' '0.1' '0.3' '1' '3' '10' '30' '100'}, ...
                  'YTickLabel', {'0.01' '0.03' '0.1' '0.3' '1' '3' '10'});
    xlabel(theAxes, 'eccentricity (degs)');
    if (displayYLabel)
        if (strcmp(theLabel, 'retinal center'))
            ylabel(theAxes, 'retinal radius (degs)');
        else
            ylabel(theAxes, 'radius (degs)');
        end
    end
    if (~isempty(theLabel))
        if (strcmp(theLabel, 'retinal center'))
            title(theAxes, 'center', 'Color', color, 'FontSize', 30);
        else
            title(theAxes, theLabel, 'Color', color, 'FontSize', 30);
        end
    end
    set(theAxes, 'XScale', 'log', 'YScale', 'log');
    grid(theAxes, 'on');
    box(theAxes, 'off');

end
