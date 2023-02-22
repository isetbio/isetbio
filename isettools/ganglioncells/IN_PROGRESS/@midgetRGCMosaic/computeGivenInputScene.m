function [responses, responseTemporalSupport, noiseFreeAbsorptionsCount] = ...
    computeGivenInputScene(obj, inputDataStruct, varargin)
% Compute the response of a midgetRGCmosaic to a scene
%
% Syntax:
%   [responses, responseTemporalSupport] = compute(obj, inputDataStruct, varargin);
%
% Description:
%    Compute the response of a midgetRGCmosaic to an input scene
%
% Inputs:
%    obj                 - A @midgetRGCMosaic object
%    inputDataStruct
%
% Outputs:
%    responses - [nTrials, nTimePoints, nConesNum] matrix
%                with noise-free midget RGC responses
%   responseTemporalSupport -

% Optional key/value pairs:
%    
%          

    % Parse optional input
    p = inputParser;
    p.addParameter('nTrials', [], @isscalar);
    p.parse(varargin{:});
    nTrials = p.Results.nTrials;

    if (isempty(nTrials))
        nTrials = 1;
    end

    % Validate the inputDataStruct
    midgetRGCMosaic.validateComputeInputDataStruct(inputDataStruct);

    % Generate the optics to use
    obj.generateOpticsAtPosition(inputDataStruct.wavefrontOpticsPositionDegs);

    % Process the null scene
    if (~isempty(inputDataStruct.theNullScene))
        % Compute the optical image of the null scene (0  contrast typically)
        obj.theCurrentOpticalImage  = oiCompute(inputDataStruct.theNullScene, obj.theCurrentOpticalImage);
        % Call the inputConeMosaic.compute() method for the current opticalImage
        noiseFreeAbsorptionsCountNull = ...
            obj.inputConeMosaic.compute(obj.theCurrentOpticalImage, 'nTrials', 1);
    end

    % Compute the optical image of the test stimulus
    obj.theCurrentOpticalImage = oiCompute(inputDataStruct.theTestScene, obj.theCurrentOpticalImage);

    % Call the inputConeMosaic.compute() method for the current opticalImage
    [inputConeMosaicSpatioTemporalActivation, ~, ~, ~, inputConeMosaicActivationTemporalSupport] = ...
        obj.inputConeMosaic.compute(obj.theCurrentOpticalImage, ...
        'nTrials', nTrials, ...
        'opticalImagePositionDegs', inputDataStruct.opticalImagePositionDegs);

    
    if (~isempty(inputDataStruct.theNullScene)) && (inputDataStruct.normalizeConeResponsesWithRespectToNullScene)
        noiseFreeAbsorptionsNormalizingResponses = 1./ noiseFreeAbsorptionsCountNull;
        idx = find(noiseFreeAbsorptionsCountNull == 0);
        noiseFreeAbsorptionsNormalizingResponses(idx) = 0;

        % Modulations from excitations
        inputConeMosaicSpatioTemporalActivation = ...
            bsxfun(@times, bsxfun(@minus, inputConeMosaicSpatioTemporalActivation, noiseFreeAbsorptionsCountNull), noiseFreeAbsorptionsNormalizingResponses);
        
    end


    % Compute the midgetRGC response by spatially pooling cone responses
    [responses, responseTemporalSupport] = obj.computeResponsesByPoolingConeResponses(...
        inputConeMosaicSpatioTemporalActivation, ...
        inputConeMosaicActivationTemporalSupport);

end
