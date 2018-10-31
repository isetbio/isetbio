function computeApertureStats(obj, plotApertureStats)
% Compute the aperture statistics in a hexagonal cone mosaic
%
% Syntax:
%    displayInfo(obj)
%
% Description:
%    Compute the statistics of the aperture diameter/area across a
%    hexagonal cone mosaic
%
% Inputs:
%    obj - The cone mosaic hex object
%    plotApertureStats - Boolean. Whether to plot a histogram of the
%    aperture sizes
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%   None.

% History:
%    10/19/2018  NPC  ISETBIO TEAM, 2015

    apertureDiametersMicrons = obj.computeApertureDiameters();
    apertureCollectingAreas = pi*(apertureDiametersMicrons/2).^2;
    obj.apertureStats.medianDiameterMicrons = median(apertureDiametersMicrons);
    obj.apertureStats.meanDiameterMicrons  = mean(apertureDiametersMicrons);
    obj.apertureStats.rangeDiameterMicrons = [min(apertureDiametersMicrons) max(apertureDiametersMicrons)];
    obj.apertureStats.medianLightCollectingArea = pi*(obj.apertureStats.medianDiameterMicrons/2)^2;
    obj.apertureStats.meanLightCollectingArea = pi*(obj.apertureStats.meanDiameterMicrons/2)^2;
    obj.apertureStats.rangeLightCollectingArea = pi*(obj.apertureStats.rangeDiameterMicrons/2).^2;
    
    if (plotApertureStats)
        hFig = figure();
        set(hFig, 'Position', [10 10 1000 440], 'Color', [1 1 1]);
    
        nBins = 15;
        subplot(1,2,1)
        histogram(apertureDiametersMicrons,nBins);
        axis 'square'
        xlabel('aperture diameter (\mum)');
        ylabel('number of cones'); 
        set(gca, 'FontSize', 16);
    
        subplot(1,2,2)
        histogram(apertureCollectingAreas, nBins);
        axis 'square'
        xlabel('Light collecting area (\mum^2)');
        ylabel('number of cones');
        set(gca, 'FontSize', 16);
    end

end