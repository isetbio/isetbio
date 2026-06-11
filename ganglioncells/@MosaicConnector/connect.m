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

    % Stage 3. Transfer input source RFs from one destinationRF (dRF1) to a
    % neighboring destination RF (dRF2), where dRF1 has at least
    % N+2 inputs, where N is the number of inputs to dRF2. This
    % transfer is done so as to minimize the combined cost for
    % dRF1+dRF2. The method called here
    % transferSourceRFsBetweenUnbalancedInputNearbyDestinationRFs()
    % is an abstract method that must be implemented by the
    % subclass so as to have specialized treatment
    obj.transferSourceRFsBetweenUnbalancedInputNearbyDestinationRFs(...
        'generateProgressVideo', false);

    % Stage 4. Transfer input source RFs from the highest input numerosity 
    % destinationRFs to any destination RFs that may have not
    % received any source RF inputs during steps1+2. 
    obj.transferSourceRFsToZeroInputDestinationRFs(...
        'generateProgressVideo', generateProgressVideo);

    % Stage 5. Swap input source RFs between nearby destination RFs 
    % of the same input numerosity so at to minimize the combined cost
    obj.swapSourceRFsBetweenNearbyDestinationRFs(...
        'generateProgressVideo', generateProgressVideo);

end % connect