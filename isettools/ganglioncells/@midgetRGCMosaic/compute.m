function [responses, responseTemporalSupport, noiseFreeAbsorptionsCount, theOpticalImage] = compute(obj, theScene, varargin)
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
    p.addParameter('withOptics', [], @isstruct);
    p.addParameter('opticalImagePositionDegs', 'mosaic-centered', @(x)(ischar(x) || (isnumeric(x)&&numel(x)==2)));
    p.addParameter('normalizeConeResponsesWithRespectToNullScene', false, @islogical);
    p.parse(varargin{:});
    

    nTrials = p.Results.nTrials;
    theNullScene = p.Results.theNullScene;
    theOpticalImage = p.Results.withOptics;
    opticalImagePositionDegs = p.Results.opticalImagePositionDegs;
    normalizeConeResponsesWithRespectToNullScene = p.Results.normalizeConeResponsesWithRespectToNullScene;

    if (isempty(nTrials))
        nTrials = 1;
    end

    if (isempty(theOpticalImage))
        % Retrieve the optics params from the RTVFTobj at the center of the mosaic
        [~, opticalPositionIndex] = min(sum((bsxfun(@minus, obj.theSamplingPositionGrid, obj.eccentricityDegs)).^2,2));
        theRTVFTobj = obj.theRetinaToVisualFieldTransformerOBJList{opticalPositionIndex};
        opticsParams = theRTVFTobj.opticsParams;
    
        fprintf('No optics were passed. Computing retinal image using optics at (%2.2f,%2.2f) degs\n', ...
            opticsParams.positionDegs(1), opticsParams.positionDegs(2));

        % Generate the OI based on the retrieved opticsParams
        oiEnsemble = obj.inputConeMosaic.oiEnsembleGenerate(opticsParams.positionDegs, ...
                        'zernikeDataBase', opticsParams.ZernikeDataBase, ...
                        'subjectID', opticsParams.testSubjectID, ...
                        'pupilDiameterMM', opticsParams.pupilDiameterMM, ...
                        'refractiveErrorDiopters', opticsParams.refractiveErrorDiopters, ...
                        'zeroCenterPSF', opticsParams.zeroCenterPSF, ...
                        'subtractCentralRefraction', opticsParams.subtractCentralRefraction, ...
                        'wavefrontSpatialSamples', opticsParams.wavefrontSpatialSamples, ...
                        'upsampleFactor', opticsParams.psfUpsampleFactor, ...
                        'warningInsteadOfErrorForBadZernikeCoeffs', true);
    
        theOpticalImage = oiEnsemble{1};
        clear 'theRTVFTobj';
    else
        fprintf('Computing retinal image using supplied optics\n');
    end

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
        obj.inputConeMosaic.compute(theOpticalImage, ...
        'nTrials', nTrials, ...
        'opticalImagePositionDegs', opticalImagePositionDegs);

    
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

    if (isempty(obj.rgcRFcenterConePoolingMatrix))
        fprintf(2,'The center and surround cone pooling matrices have not yet been set (no center/surround RF and no overlap) !!\n');
        fprintf(2,'Using the rgcRFcenterConeConnectivityMatrix (zero overlap) instead to compute RF center responses only !!\n');

        % Compute the response of each mRGC
        parfor iRGC = 1:mRGCsNum
            % Retrieve the center cone indices & weights
            centerConnectivityVector = full(squeeze(obj.rgcRFcenterConeConnectivityMatrix(:, iRGC)));
            centerConeIndices = find(centerConnectivityVector > 0.0001);
            centerConeWeights = reshape(centerConnectivityVector(centerConeIndices), [1 1 numel(centerConeIndices)]);
    
            % Sum the weighted cone responses pooled by the RF center only
            centerActivation = sum(bsxfun(@times, noiseFreeAbsorptionsCount(:,:,centerConeIndices), centerConeWeights),3);
            responses(:,:,iRGC) = centerActivation;
        end

    else

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

end
