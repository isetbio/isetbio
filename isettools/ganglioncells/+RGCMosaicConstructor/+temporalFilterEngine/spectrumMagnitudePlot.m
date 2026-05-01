%
% RGCMosaicConstructor.temporalFilterEngine.spectrumMagnitudePlot
%

function pHandle = spectrumMagnitudePlot(ax, temporalFrequencySupportHz, theTTF, theMarker, ...
    noXLabel, noYLabel, theColor, theLabel)

    theMagnitudeSpectrum = abs(theTTF);
    theMagnitudeSpectrum = theMagnitudeSpectrum / max(theMagnitudeSpectrum(:));

    if (strcmp(theMarker, '-'))
         pHandle = plot(ax, temporalFrequencySupportHz, theMagnitudeSpectrum, '-', ...
             'Color', theColor, 'LineWidth', 1.5);
    else
        pHandle = scatter(ax, temporalFrequencySupportHz, theMagnitudeSpectrum, 144, ...
             'Marker', theMarker, 'Color', theColor, 'LineWidth', 1.5, 'MarkerEdgeColor', theColor*0.5, ...
             'MarkerFaceColor', theColor, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
    end

    set(ax, 'XScale', 'log', 'XLim', [0.3 200], 'XTick', [0.1 0.3 1 3 10 30 100]);
    set(ax, 'YScale', 'log', 'YLim', [0.001 1.0], 'YTick', [0.001 0.01 0.1 1], 'YTickLabel', {'.001' '.01', '.1', '1'})
    set(ax, 'FontSize', 16);
    
    if (~noXLabel)
        xlabel(ax, 'frequency (Hz)')
    end

     if (~noYLabel)
        ylabel(ax, 'magnitude (normalized)')
     end

    if (~isempty(theLabel))
        title(ax, theLabel);
    end
    
end
