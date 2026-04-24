%
% RGCMosaicConstructor.temporalFilterEngine.spectrumPhasePlot
%

function spectrumPhasePlot(ax, temporalFrequencySupportHz, theTTF, noXLabel, noYLabel, theColor, theLabel)

    plot(ax, temporalFrequencySupportHz, unwrap(angle(theTTF))/pi*180, 'ko', 'LineWidth', 1.0, 'MarkerEdgeColor', theColor, 'MarkerFaceColor', 0.5*theColor + 0.5*[1 1 1]);
    set(ax, 'XScale', 'log', 'XLim', [0.3 300], 'XTick', [0.1 0.3 1 3 10 30 100]);

    set(ax, 'FontSize', 16);
    if (~noXLabel)
        xlabel(ax, 'frequency (Hz)')
    end

     if (~noYLabel)
        ylabel(ax, 'phase (degs)')
     end
    
    if (~isempty(theLabel))
        title(ax, theLabel);
    end

end
