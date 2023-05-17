function outContainer = combineContainersMat(theCellArrayOfContainers)
% Consolidate a cell array of containers of matrices into a container of a matrix.
% 
% Syntax:
%    outContainer = combineContainersMat(theCellArrayOfContainers)
%
% Description:
%    Take a cell array of containers, each of which contains a matrix of
%    the same dimension, and and returns a single ontainer whose contents
%    are a matrix that has the individual matrices packed into it.
%
%    The resulting matrix has one more dimension than the individual
%    matrices, with the first dimension indexing which of the passed cell
%    arrays the matrix in the rest of the dimensions came from.
%    
%    This operates over all the flags of the passed containers, and check
%    that the flags are matched across the containers passed.
%
%    The matrices have to be the same size within flag, but can differ
%    across flags.
%
% Inputs:
%    theCellArrayOfContainers      - The passed cell array of containers
%
% Outputs:
%    outContainer                  - The consolidated container of cell
%                                    arrays.
%
% Optional key/value pairs:
%    None
%
% See also: combineContainers.
%

% History:
%    12/13/21  dhb  Wrote it.
%    31/03/23  NPC  Update it to properly deal with containers containing
%                   multiple trials

    % Ensure that all containers have the same flags
    validateContainerKeys(theCellArrayOfContainers);
    
    % Initialize output container
    outContainer = containers.Map();
    
    % Fill the entries of the outContainer
    theKeys = keys(theCellArrayOfContainers{1});
    for iEntry = 1:numel(theKeys)
        theKey = theKeys{iEntry};
        theKeyValuesAcrossAllContainers = cell(1,numel(theCellArrayOfContainers));
        for iContainer = 1:numel(theCellArrayOfContainers)
            % Allocate matrix for this key
            if (iContainer == 1)
                theSpatioTemporalResponses = theCellArrayOfContainers{iContainer}(theKey);
                assert(ndims(theSpatioTemporalResponses)==3, 'The spatiotemporal responses matrix must have 3 dimensions');

                theMatSize = size(theSpatioTemporalResponses);
                nTrials = theMatSize(1);
                theSingleTrialSpatioTemporalResponseSize = prod(theMatSize(2:3));

                theMat = zeros([numel(theCellArrayOfContainers)*nTrials theSingleTrialSpatioTemporalResponseSize]);
            else
                if (any(size(theCellArrayOfContainers{iContainer}(theKey)) ~= theMatSize))
                    error('Mismatch of matrix sizes within a key');
                end
            end
            % Retrieve the spatiotemporal responses
            theSpatioTemporalResponses = theCellArrayOfContainers{iContainer}(theKey);

            % Reshape into a row vector
            for iTrial = 1:nTrials
                theSingleTrialSpatioTemporalResponse = theSpatioTemporalResponses(iTrial,:,:);
                theMat((iContainer-1)*nTrials+iTrial,:) = reshape(theSingleTrialSpatioTemporalResponse, [1 theSingleTrialSpatioTemporalResponseSize]);
            end

        end
        outContainer(theKey) = theMat;
    end
end

function validateContainerKeys(theCellArrayOfContainers)
% Ensure that all containers have the same flags

    for iContainer = 1:numel(theCellArrayOfContainers)
        theContainerKeys = keys(theCellArrayOfContainers{iContainer});
        if (iContainer == 1)
            theCommonKeys = theContainerKeys;
        else
            % Check that all keys in this container are also in theCommonKeys
            for iKey = 1:numel(theContainerKeys)
                assert(ismember(theContainerKeys{iKey}, theCommonKeys), 'key ''%s'' in the %d container does not exist in the first container', ...
                    theContainerKeys{iKey}, iContainer);
            end
            % Check that all keys in theCommonKeys are also keys in this container
            for iKey = 1:numel(theCommonKeys)
                assert(ismember(theCommonKeys{iKey}, theContainerKeys), 'key ''%s'' in the first container does not exist in the %d container', ...
                    theCommonKeys{iKey}, iContainer);
            end
        end
    end
end
