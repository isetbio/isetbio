function [responses, responseTemporalSupport] = ...
    compute(obj, inputDataStruct, varargin)
% Compute the response of a midgetRGCmosaic to a scene
%
% Syntax:
%   [responses, responseTemporalSupport] = compute(obj, inputDataStruct, varargin);
%
% Description:
%    Compute the response of a midgetRGCmosaic to an input. This is a
%    gateway function that diverts to one of the more specialized compute
%    methods, depending on the type of inputData
%
% Inputs:
%    obj                 - A @midgetRGCMosaic object
%    inputDataStruct     -
%
% Outputs:
%    responses - [nTrials, nTimePoints, nRGCsNum] matrix
%                with noise-free midget RGC responses
%   responseTemporalSupport -

% Optional key/value pairs:
%    
%          

    % Validate inputDataStruct
    midgetRGCMosaic.validateComputeInputDataStruct(inputDataStruct);

    % Parse optional input
    p = inputParser;
    p.addParameter('nTrials', [], @isscalar);
    p.parse(varargin{:});
    nTrials = p.Results.nTrials;

    switch (inputDataStruct.type)
        case midgetRGCMosaic.SCENE_COMPUTE_INPUT_DATA_TYPE 
            [responses, responseTemporalSupport] = ...
                obj.computeGivenInputScene(inputDataStruct, ...
                    'nTrials', nTrials);

        case midgetRGCMosaic.CONE_MOSAIC_RESPONSE_COMPUTE_INPUT_DATA_TYPE 
            [responses, responseTemporalSupport] = ...
                obj.computeGivenInputConeMosaicActivation(inputDataStruct, ...
                    'nTrials', nTrials);

    end

end
