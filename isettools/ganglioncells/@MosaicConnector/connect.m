function connect(obj, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('generateProgressVideo', false, @islogical);
    p.parse(varargin{:});
    generateProgressVideo = p.Results.generateProgressVideo;


    % Empty the intermediateFigureHandles cell array
    obj.intermediateFigureHandles = {};
    % Empty the intermediate meta data struct cell array
    obj.intermediateMetaDataStructs = {};


    % Stage 1. Connect mosaics based on the sourceToDestinationDensityRatio
    obj.connectSourceRFsToDestinationRFsBasedOnLocalDensities();
 
    % Stage 2. Connect unconnected sourceRFs to their closest destination RF
    obj.connectUnconnectedSourceRFsToClosestDestinationRF();

    % Stage 3. Transfer input source RFs from the highest input numerosity 
    % destinationRFs to any destination RFs with ZERO inputs 
    obj.transferSourceRFsToZeroInputDestinationRFs(...
        'generateProgressVideo', generateProgressVideo);

    % Stage 4. Transfer input source RFs from the highest input numerosity 
    % destinationRFs to any nearby destinationRFs with LESS # of inputs
    obj.transferSourceRFsBetweenNearbyDestinationRFs(...
        'generateProgressVideo', generateProgressVideo, ...
        'swapSourceRFsInsteadOffUnilateralTransfer', false);

    % Stage 5. Swap input source RFs between nearby destination RFs 
    % of the same input numerosity so at to minimize the combined cost
    obj.transferSourceRFsBetweenNearbyDestinationRFs(...
        'generateProgressVideo', generateProgressVideo, ...
        'swapSourceRFsInsteadOffUnilateralTransfer', true);

    
    if (obj.wiringParams.spatialChromaticUniformityTradeoff == 1)
        % Repeat Stage 4 is we are at max spatial compactness.
        % Some times this helps remove minor spatial compactness anomalies introduced at Stage 5
        obj.transferSourceRFsBetweenNearbyDestinationRFs(...
            'generateProgressVideo', generateProgressVideo, ...
            'swapSourceRFsInsteadOffUnilateralTransfer', false);
    end

end % connect