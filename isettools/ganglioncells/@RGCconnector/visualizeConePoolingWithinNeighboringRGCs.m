function [hFig,visualizedNeighborsNum] = visualizeConePoolingWithinNeighboringRGCs(obj, iRGC, varargin)
    
     % Parse input
    p = inputParser;
    p.addParameter('visualizedFieldOfViewMicrons', [], @(x)((isempty(x))||(isscalar(x))));
    p.addParameter('visualizedConesNum', [], @(x)((isempty(x))||(isscalar(x))));

    p.parse(varargin{:});
    visualizedFieldOfViewMicrons = p.Results.visualizedFieldOfViewMicrons;
    visualizedConesNum = p.Results.visualizedConesNum;

    
    % Find the closest neighboring RGC
    nearbyRGCindices = obj.neihboringRGCindices(iRGC);
    visualizedNeighborsNum = numel(nearbyRGCindices);
    if (isempty(nearbyRGCindices))
        return;
    end
    
    
    
    % Set-up figure
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'colsNum', max([6 numel(nearbyRGCindices)]), ...
       'rowsNum', 3, ...
       'heightMargin',  0.09, ...
       'widthMargin',    0.04, ...
       'leftMargin',     0.03, ...
       'rightMargin',    0.00, ...
       'bottomMargin',   0.05, ...
       'topMargin',      0.02);
   
   
    hFig = figure(); clf;
    set(hFig, 'Position', [10 10 2300 1000], 'Color', [1 1 1]);
    
    % Plot the main RGC wiring
    for iCol = 1:numel(nearbyRGCindices)
        ax = subplot('Position', subplotPosVectors(1,iCol).v);
        [xSupport, ySupport, rfProfile2DmainRGC, xTicks, yTicks] = ...
            obj.visualizeConeAperturePooling(iRGC, ax, [], ...
            [], [], [], [], visualizedFieldOfViewMicrons, visualizedConesNum, 'black');
        
        % Plot the nearby RGC wiring
        ax = subplot('Position', subplotPosVectors(2,iCol).v);
        [~, ~, rfProfile2DnearbyRGC, ~, ~] = ...
            obj.visualizeConeAperturePooling(nearbyRGCindices(iCol), ax, [],  ...
            xSupport, ySupport, xTicks, yTicks, visualizedFieldOfViewMicrons, visualizedConesNum, 'red');
        
        ax = subplot('Position', subplotPosVectors(3,iCol).v);
        obj.visualizeRFoverlap2D(iRGC, nearbyRGCindices(iCol), rfProfile2DmainRGC, rfProfile2DnearbyRGC, ... 
            xSupport, ySupport, xTicks, yTicks, ax);
    end
    
    
end

