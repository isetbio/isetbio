function adjustRFoverlap(obj, overlapRatio)

    hFig = figure(6969); clf;
    ax = subplot(2,2,1);
    [~,~,XLims, YLims] = obj.theMosaicConnectorOBJ.visualizeInputLattices(...
        'figureHandle', hFig, ...
        'axesHandle', ax, ...
        'thetaSamples', 30);
    set(ax, 'FontSize', 16);

    ax = subplot(2,2,2);
    obj.theMosaicConnectorOBJ.visualizeDestinationLatticePooling(...
        'figureHandle', hFig, ...
        'axesHandle', ax, ...
        'titleString', 'BeforeDiverge', ...
        'titleWithPoolingStats', true, ...
        'XLims', XLims, ...
        'YLims', YLims);

    % Update the rfOverlapRatio property
    obj.rfOverlapRatio = overlapRatio;

    % Diverge input cones to multiple nearby midget RGCs (RF overlap)
    obj.theMosaicConnectorOBJ.divergeSourceRFsToNearbyDestinationRFs(...
        'destinationRFoverlapRatio', obj.rfOverlapRatio);

    % Set the rf center connectivity matrix
    obj.rgcRFcenterConeConnectivityMatrix = obj.theMosaicConnectorOBJ.connectivityMatrix;
                
    % Set the rf positions, in microns
    obj.rgcRFpositionsMicrons = obj.theMosaicConnectorOBJ.destinationRFcentroidsFromInputs;

    % Set the rf spacings, in microns
    obj.rgcRFspacingsMicrons = obj.theMosaicConnectorOBJ.destinationRFspacingsFromCentroids;
   
    % Set the rf positions in degrees
    obj.rgcRFpositionsDegs = obj.inputConeMosaic.customMMsToDegsConversionFunction(obj.rgcRFpositionsMicrons/1000);
    
    % Set the rf spacings in degrees
    eccentricityMicrons = sqrt(sum((obj.rgcRFpositionsMicrons).^2,2));
    obj.rgcRFspacingsDegs = obj.inputConeMosaic.sizeMicronsToSizeDegreesForCmosaic(obj.rgcRFspacingsMicrons, eccentricityMicrons');

    % Crop RGCs in the periphery
    obj.cropRGCsOnTheBorder();

    ax = subplot(2,2,4);
    obj.theMosaicConnectorOBJ.visualizeDestinationLatticePooling(...
        'figureHandle', hFig, ...
        'axesHandle', ax, ...
        'titleString', 'AfterDiverge', ...
        'titleWithPoolingStats', true, ...
        'XLims', XLims, ...
        'YLims', YLims);
    
end


