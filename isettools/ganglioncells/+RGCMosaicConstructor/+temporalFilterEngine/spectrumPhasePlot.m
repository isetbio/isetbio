%
% RGCMosaicConstructor.temporalFilterEngine.spectrumPhasePlot
%

function pHandle = spectrumPhasePlot(ax, temporalFrequencySupportHz, theTTF, theMarker, noXLabel, noYLabel, theColor, theLabel)

    theUnwrappedPhaseRadians = unwrap(angle(theTTF));

    while (theUnwrappedPhaseRadians(1) > 0.5*pi)
        theUnwrappedPhaseRadians = theUnwrappedPhaseRadians - 2*pi;
    end

    if (strcmp(theMarker, '-'))
         pHandle = plot(ax, temporalFrequencySupportHz, theUnwrappedPhaseRadians/pi, '-', ...
             'Color', theColor, 'LineWidth', 1.5);
    else
        pHandle = scatter(ax, temporalFrequencySupportHz, theUnwrappedPhaseRadians/pi, 144, ...
             'Marker', theMarker, 'Color', theColor, 'LineWidth', 1.5, 'MarkerEdgeColor', theColor*0.5, ...
             'MarkerFaceColor', theColor, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
    end


    set(ax, 'XScale', 'log', 'XLim', [0.3 200], 'XTick', [0.1 0.3 1 3 10 30 100], 'YLim', [-7 0.5], 'YTick', -10:1:10);

    set(ax, 'FontSize', 16);
    if (~noXLabel)
        xlabel(ax, 'frequency (Hz)')
    end

     if (~noYLabel)
        ylabel(ax, 'phase (pi radians)')
     end
    
    if (~isempty(theLabel))
        title(ax, theLabel);
    end

end
