% Method to return the cone-to-RGC density ratio for all RGCs
function densityRatios = coneToRGCDensityRatios(obj)

    xyRange = obj.maxRFpositionMicrons - obj.minRFpositionMicrons;
    N = 32;
    if (xyRange(1) > xyRange(2))
        nY = N;
        nX = round(nY * xyRange(1)/xyRange(2));
    else
        nX = N;
        nY = round(nX * xyRange(2)/xyRange(1));
    end
    
    % Sampling grid
    x = linspace(obj.minRFpositionMicrons(1), obj.maxRFpositionMicrons(1), nX);
    y = linspace(obj.minRFpositionMicrons(2), obj.maxRFpositionMicrons(2), nY);
    gridPositions{1} = x(:); 
    gridPositions{2} = y(:);
    
    % Compute cone mosaic density map
    coneDensity2DMap = cMosaic.densityMap(...
        obj.inputConeMosaic.coneRFpositionsMicrons, ...
        obj.inputConeMosaic.coneRFspacingsMicrons, gridPositions);
    
    % Compute rgcmosaic density map
    rgcDensity2DMap = cMosaic.densityMap(...
        obj.rgcRFpositionsMicrons, obj.rgcRFspacingsMicrons, gridPositions);
    
    % Ratios
    coneToRGCratiosGrid = (coneDensity2DMap ./ rgcDensity2DMap);
    
    % Compute gridded interpolant for the original cone positions
    interpolationMethod = 'linear';
    extrapolationMethod = 'nearest';
    F = griddedInterpolant(gridPositions, coneToRGCratiosGrid', interpolationMethod, extrapolationMethod);
    
    % Compute ratio for each RGC
    densityRatios = F(obj.rgcRFpositionsMicrons);
    
end

