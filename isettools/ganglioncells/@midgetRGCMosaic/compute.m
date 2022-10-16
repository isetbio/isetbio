function [responses, responseTemporalSupport] = compute(obj, theScene, varargin)
% Compute the response of a midgetRGCmosaic to a scene
%
% Syntax:
%   [responses, responseTemporalSupport] = compute(obj, theScene);
%
% Description:
%    Compute the response of a midgetRGCmosaic to an optical image
%
% Inputs:
%    obj                 - A @midgetRGCMosaic object
%    theScene            - A scene
%
% Outputs:
%    responses - [nTrials, nTimePoints, nConesNum] matrix
%                with noise-free midget RGC responses
%   responseTemporalSupport -

% Optional key/value pairs:
%    
%          

    % Parse input
    p = inputParser;
    p.addParameter('nTrials', [], @isscalar);
    p.addParameter('theNullScene', [], @isstruct);
    p.addParameter('normalizeConeResponsesWithRespectToNullScene', false, @islogical);
    p.parse(varargin{:});
    

    nTrials = p.Results.nTrials;
    theNullScene = p.Results.theNullScene;
    normalizeConeResponsesWithRespectToNullScene = p.Results.normalizeConeResponsesWithRespectToNullScene;

    if (isempty(nTrials))
        nTrials = 1;
    end

    % Retrieve the optics
    % Note: if the obj.theOpticsPositionGrid contains more than one
    % position, we will have to generate multiple optical images
    % more than 1 
    opticalPositionsNum = size(obj.theOpticsPositionGrid,1);
    opticalPositionIndex = 1;
    if (opticalPositionsNum > 1)
        fprintf(2,'Computing with multiple optical images is not yet supported. Using the first one.\n');
    end

    % Retrieve the optics from the first RTVFTobj
    theRTVFTobj = obj.theRetinaToVisualFieldTransformerOBJList{opticalPositionIndex};
    theOpticalImage = theRTVFTobj.theOI;
    clear 'theRTVFTobj';

    % Process the null scene
    if (~isempty(theNullScene))
        % Compute the optical image of the null scene (0  contrast typically)
        theOpticalImage = oiCompute(theNullScene, theOpticalImage);
        % Call the inputConeMosaic.compute() method for the current opticalImage
        [noiseFreeAbsorptionsCountNull, ~, ...
         photoCurrentsNull, ~, ~] = ...
            obj.inputConeMosaic.compute(theOpticalImage, 'nTrials', 1);
    end

    % Compute the optical image of the test stimulus
    theOpticalImage = oiCompute(theScene, theOpticalImage);

    % Call the inputConeMosaic.compute() method for the current opticalImage
    [noiseFreeAbsorptionsCount, noisyAbsorptionsCountInstances, ...
     photoCurrents, photoCurrentInstances, responseTemporalSupport] = ...
        obj.inputConeMosaic.compute(theOpticalImage,'nTrials', nTrials);

    
    if (~isempty(theNullScene)) && (normalizeConeResponsesWithRespectToNullScene)
        noiseFreeAbsorptionsCount = ...
            (bsxfun(@minus, noiseFreeAbsorptionsCount, noiseFreeAbsorptionsCountNull)) ./ noiseFreeAbsorptionsCountNull;
        noisyAbsorptionsCountInstances = ...
            (bsxfun(@minus, noisyAbsorptionsCountInstances, noiseFreeAbsorptionsCountNull)) ./ noiseFreeAbsorptionsCountNull;
    end

    nTrials = size(noiseFreeAbsorptionsCount,1);
    nTimePoints = size(noiseFreeAbsorptionsCount,2);
    mRGCsNum = size(obj.rgcRFcenterConeConnectivityMatrix,2);
    
    % Allocate memory for the responses
    responses = zeros(nTrials, nTimePoints, mRGCsNum);

    % Compute the response of each mRGC
    parfor iRGC = 1:mRGCsNum
        % Retrieve the center cone indices & weights
        centerConnectivityVector = full(squeeze(obj.rgcRFcenterConePoolingMatrix(:, iRGC)));
        centerConeIndices = find(centerConnectivityVector > 0.0001);
        centerConeWeights = reshape(centerConnectivityVector(centerConeIndices), [1 1 numel(centerConeIndices)]);

        % Retrieve the surround cone indices & weights
        surroundConnectivityVector = full(squeeze(obj.rgcRFsurroundConePoolingMatrix(:, iRGC)));
        surroundConeIndices = find(surroundConnectivityVector > 0.0001);
        surroundConeWeights = reshape(surroundConnectivityVector(surroundConeIndices), [1 1 numel(surroundConeIndices)]);

        % Sum the weighted cone responses pooled by the RF center and RF surround
        centerActivation = sum(bsxfun(@times, noiseFreeAbsorptionsCount(:,:,centerConeIndices), centerConeWeights),3);
        surroundActivation = sum(bsxfun(@times, noiseFreeAbsorptionsCount(:,:,surroundConeIndices), surroundConeWeights),3);

        % Composite response: centerActivation - surroundActivation
        responses(:,:,iRGC) = centerActivation - surroundActivation;
    end

end
