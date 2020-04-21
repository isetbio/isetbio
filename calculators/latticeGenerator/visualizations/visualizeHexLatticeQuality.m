function visualizeHexLatticeQuality(histogramData, minQualityValue)
        
    hFig = figure(2);
    bar(histogramData.x, histogramData.y, 'FaceColor', [0.5 0.5 0.5]); hold on;
    plot(minQualityValue*[1 1], [0 max(histogramData.y)], 'r');
    grid on
    xlabel('hex-index $\left(\displaystyle 2 r_{ins} / r_{cir} \right)$', 'Interpreter', 'latex');
  
end

