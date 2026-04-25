%
% RGCMosaicConstructor.temporalFilterEngine.spectrumMagnitudePlot
%

function spectrumMagnitudePlot(ax, temporalFrequencySupportHz, theTTF, theMarker, ...
    noXLabel, noYLabel, theColor, theLabel)

    plot(ax, temporalFrequencySupportHz, abs(theTTF), theMarker, ...
         'Color', theColor, 'LineWidth', 1.5, 'MarkerEdgeColor', theColor, 'MarkerFaceColor', 0.5*theColor + 0.5*[1 1 1]);
    set(ax, 'XScale', 'log', 'XLim', [0.3 300], 'XTick', [0.1 0.3 1 3 10 30 100])
    set(ax, 'FontSize', 16);
    
    if (~noXLabel)
        xlabel(ax, 'frequency (Hz)')
    end

     if (~noYLabel)
        ylabel(ax, '|H|')
     end

    if (~isempty(theLabel))
        title(ax, theLabel);
    end
    
end
