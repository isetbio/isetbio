function [responses, responseTemporalSupport] = ...
    computeGivenInputConeMosaicActivation(obj, inputDataStruct, varargin)
% Compute the response of a midgetRGCmosaic to a scene
%
% Syntax:
%   [responses, responseTemporalSupport] = compute(obj, varargin);
%
% Description:
%    Compute the response of a midgetRGCmosaic to an input scene
%
% Inputs:
%    obj                 - A @midgetRGCMosaic object
%    inputDataStruct
%
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


    % Compute the midgetRGC response by spatially pooling cone responses
    [responses, responseTemporalSupport] = obj.computeResponsesByPoolingConeResponses(...
        inputDataStruct.inputConeMosaicSpatioTemporalActivation, ...
        inputDataStruct.inputConeMosaicActivationTemporalSupport ...
        );

end
