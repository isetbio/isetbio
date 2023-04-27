function cropToSizeAtEccentricity(obj, sizeDegs, eccentricityDegs, varargin)
    % Parse input
    p = inputParser;
    p.addRequired('sizeDegs', @isnumeric);
    p.addRequired('eccentricityDegs', @isnumeric);
    p.addParameter('name', '', @ischar);
    p.addParameter('visualizeSpatialRelationshipToSourceMosaic', false, @islogical);
    p.addParameter('beVerbose', false, @islogical);

    % Parse input
    p.parse(sizeDegs, eccentricityDegs, varargin{:});

    assert(~isempty(obj.rgcRFcenterConePoolingMatrix), 'This is not a compute-ready mosaic.');
    assert(~isempty(obj.rgcRFsurroundConePoolingMatrix), 'This is not a compute-ready mosaic.');
    assert(numel(eccentricityDegs)==2, 'eccentricity must be a 2-element vector');
    assert(numel(sizeDegs)==2, 'size must be a 2-element vector');


    name = p.Results.name;
    visualizeSpatialRelationshipToSourceMosaic = p.Results.visualizeSpatialRelationshipToSourceMosaic;
    beVerbose = p.Results.beVerbose;


    % Define theROI based on the passed eccentricity and size
    theROI = regionOfInterest(...
        'geometryStruct', struct(...
            'units', 'degs', ...
            'shape', 'rect', ...
            'center', eccentricityDegs, ...
            'width', sizeDegs(1), ...
            'height', sizeDegs(2), ...
            'rotation', 0.0...
        ));

    % Find the indices of the RGCs whose position is within theROI
    keptRGCindices = theROI.indicesOfPointsInside(obj.rgcRFpositionsDegs);

    % Find cone indices to keep
    keptConeIndices = [];
    for iRGC = 1:numel(keptRGCindices)
        coneIndices = find(squeeze(obj.rgcRFsurroundConePoolingMatrix(:, keptRGCindices(iRGC)))>0.000001);
        keptConeIndices = unique(cat(1,keptConeIndices(:),coneIndices(:)));
    end

    % Determine the indices of cones to keep
    xPos = obj.inputConeMosaic.coneRFpositionsDegs(keptConeIndices,1);
    yPos = obj.inputConeMosaic.coneRFpositionsDegs(keptConeIndices,2);
    xMinPos = min(xPos(:));
    xMaxPos = max(xPos(:));
    yMinPos = min(yPos(:));
    yMaxPos = max(yPos(:));

    keptAllConeTypeIndices = find(...
        (obj.inputConeMosaic.coneRFpositionsDegs(:,1) >= xMinPos) & ...
        (obj.inputConeMosaic.coneRFpositionsDegs(:,1) <= xMaxPos) & ...
        (obj.inputConeMosaic.coneRFpositionsDegs(:,2) >= yMinPos) & ...
        (obj.inputConeMosaic.coneRFpositionsDegs(:,2) <= yMaxPos) ...
        );

    % Update obj.inputConeMosaic to only include the keptConeIndices
    obj.inputConeMosaic.cropMosaicToIncluceConesWithIndices(keptAllConeTypeIndices);

    % Update the rgcRFpositions and spacings of the cropped mRGCMosaic 
    obj.rgcRFpositionsMicrons = obj.rgcRFpositionsMicrons(keptRGCindices,:);
    obj.rgcRFpositionsDegs = obj.rgcRFpositionsDegs(keptRGCindices,:);
    obj.rgcRFspacingsMicrons = obj.rgcRFspacingsMicrons(1,keptRGCindices);
    obj.rgcRFspacingsDegs = obj.rgcRFspacingsDegs(1,keptRGCindices);

    % Crop obj.rgcRFcenterConePoolingMatrix to generate croppedCenterConePoolingMatrix
    croppedCenterConePoolingMatrix = obj.rgcRFcenterConePoolingMatrix(keptAllConeTypeIndices, keptRGCindices);

    % Crop obj.rgcRFsurroundConePoolingMatrix to generate croppedSurroundConePoolingMatrix
    croppedSurroundConePoolingMatrix = obj.rgcRFsurroundConePoolingMatrix(keptAllConeTypeIndices, keptRGCindices);

    % Update cone pooling matrices with the cropped versions
    obj.rgcRFcenterConePoolingMatrix = croppedCenterConePoolingMatrix;
    obj.rgcRFsurroundConePoolingMatrix = croppedSurroundConePoolingMatrix;

    % Update the eccentricity and size of the cropped mRGCMosaic 
    obj.eccentricityDegs = eccentricityDegs;
    obj.eccentricityMicrons = obj.inputConeMosaic.distanceDegreesToDistanceMicronsForCmosaic(eccentricityDegs);
    obj.sizeDegs = sizeDegs;

    % Reset the visualizationCache
    obj.visualizationCache = [];
    
end
