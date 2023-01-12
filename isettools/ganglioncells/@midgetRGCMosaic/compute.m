function [responses, responseTemporalSupport, noiseFreeAbsorptionsCount] = compute(obj, varargin)
% Compute the response of a midgetRGCmosaic to a scene
%
% Syntax:
%   [responses, responseTemporalSupport] = compute(obj, varargin);
%
% Description:
%    Compute the response of a midgetRGCmosaic to an optical image
%
% Inputs:
%    obj                 - A @midgetRGCMosaic object
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

    % The following parameters only make sense when we are asked to compute
    % midgetRGCMosaic responses to scenes. Probably should combine all
    % these in a struct
    p.addParameter('theTestScene', [], @isstruct);
    p.addParameter('theNullScene', [], @isstruct);
    p.addParameter('withWavefronOpticsAtPositionDegs', obj.eccentricityDegs, @(x)((isnumeric(x))&&(numel(x)==2)))
    p.addParameter('opticalImagePositionDegs', 'mosaic-centered', @(x)(ischar(x) || (isnumeric(x)&&numel(x)==2)));
    p.addParameter('normalizeConeResponsesWithRespectToNullScene', false, @islogical);

    % The following parameters only make sense when we are asked to compute
    % midgetRGCMosaic responses from cone mosaic responses which were precomputed 
    % by the user (using the obj.inputConeMosaic)
    


    p.parse(varargin{:});
    

    nTrials = p.Results.nTrials;
    theScene = p.Results.theTestScene;
    theNullScene = p.Results.theNullScene;
    wavefrontOpticsPositionDegs = p.Results.withWavefronOpticsAtPositionDegs;
    opticalImagePositionDegs = p.Results.opticalImagePositionDegs;
    normalizeConeResponsesWithRespectToNullScene = p.Results.normalizeConeResponsesWithRespectToNullScene;

    if (isempty(nTrials))
        nTrials = 1;
    end

    % Generate the optics to use
    obj.generateOpticsAtPosition(wavefrontOpticsPositionDegs);

    % Process the null scene
    if (~isempty(theNullScene))
        % Compute the optical image of the null scene (0  contrast typically)
        obj.theCurrentOpticalImage  = oiCompute(theNullScene, obj.theCurrentOpticalImage);
        % Call the inputConeMosaic.compute() method for the current opticalImage
        [noiseFreeAbsorptionsCountNull, ~, ...
         photoCurrentsNull, ~, ~] = ...
            obj.inputConeMosaic.compute(obj.theCurrentOpticalImage, 'nTrials', 1);
    end

    % Compute the optical image of the test stimulus
    obj.theCurrentOpticalImage = oiCompute(theScene, obj.theCurrentOpticalImage);

    % Call the inputConeMosaic.compute() method for the current opticalImage
    [noiseFreeAbsorptionsCount, noisyAbsorptionsCountInstances, ...
     photoCurrents, photoCurrentInstances, responseTemporalSupport] = ...
        obj.inputConeMosaic.compute(obj.theCurrentOpticalImage, ...
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

    if (isempty(obj.rgcRFcenterConePoolingMatrix))
        mRGCsNum = size(obj.rgcRFcenterConeConnectivityMatrix,2);
    else
        mRGCsNum = size(obj.rgcRFcenterConePoolingMatrix,2);
    end

    % Allocate memory for the responses
    responses = zeros(nTrials, nTimePoints, mRGCsNum);

    if (isempty(obj.rgcRFcenterConePoolingMatrix))
        fprintf(2,'The center and surround cone pooling matrices have not yet been set (no center/surround RF and no overlap) !!\n');
        fprintf(2,'Using the rgcRFcenterConeConnectivityMatrix instead to compute RF center responses only !!\n');

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
