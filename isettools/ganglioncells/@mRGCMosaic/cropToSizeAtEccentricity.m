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

    % Find the RGCs whose position is within theROI
    keptRGCindices = theROI.indicesOfPointsInside(obj.rgcRFpositionsDegs);

    % Update the eccentricity and size of the cropped mRGCMosaic 
    obj.eccentricityDegs = eccentricityDegs;
    obj.eccentricityMicrons = obj.inputConeMosaic.distanceDegreesToDistanceMicronsForCmosaic(eccentricityDegs);
    obj.sizeDegs = sizeDegs;

    % Update the inputConeMosaic of the cropped mRGCMosaic 
    % For now we are keeping the entire inputConeMosaic
    keptConeIndices = 1:obj.inputConeMosaic.conesNum;

    % Update the rgcRFpositions and spacings of the cropped mRGCMosaic 
    obj.rgcRFpositionsMicrons = obj.rgcRFpositionsMicrons(keptRGCindices,:);
    obj.rgcRFpositionsDegs = obj.rgcRFpositionsDegs(keptRGCindices,:);
    obj.rgcRFspacingsMicrons = obj.rgcRFspacingsMicrons(1,keptRGCindices);
    obj.rgcRFspacingsDegs = obj.rgcRFspacingsDegs(1,keptRGCindices);

    % Crop obj.rgcRFcenterConePoolingMatrix to generate croppedCenterConePoolingMatrix
    croppedCenterConePoolingMatrix = obj.rgcRFcenterConePoolingMatrix(keptConeIndices, keptRGCindices);

    % Crop obj.rgcRFsurroundConePoolingMatrix to generate croppedSurroundConePoolingMatrix
    croppedSurroundConePoolingMatrix = obj.rgcRFsurroundConePoolingMatrix(keptConeIndices, keptRGCindices);

    % Update cone pooling matrices with the cropped versions
    obj.rgcRFcenterConePoolingMatrix = croppedCenterConePoolingMatrix;
    obj.rgcRFsurroundConePoolingMatrix = croppedSurroundConePoolingMatrix;

    % Reset the visualizationCache
    obj.visualizationCache = [];
    
end
