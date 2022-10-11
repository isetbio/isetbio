function [responses, responseTemporalSupport] = compute(obj, oi, varargin)
% Compute the response of a midgetRGCmosaic to an optical image
%
% Syntax:
%   [responses, responseTemporalSupport] = compute(obj, oi);
%
% Description:
%    Compute the response of a midgetRGCmosaic to an optical image
%
% Inputs:
%    obj                 - A @midgetRGCMosaic object
%    oi                  - An optical image, or oiSequence. 
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
    p.parse(varargin{:});
    

    nTrials = p.Results.nTrials;
    if (isempty(nTrials))
        nTrials = 1;
    end

    % Call the inputConeMosaic.comptute() method
    [noiseFreeAbsorptionsCount, noisyAbsorptionInstances, ...
     photoCurrents, photoCurrentInstances, responseTemporalSupport] = obj.inputConeMosaic.compute(oi, ...
        'nTrials', nTrials);

    nTrials = size(noiseFreeAbsorptionsCount,1);
    nTimePoints = size(noiseFreeAbsorptionsCount,2);
    conesNum = size(noiseFreeAbsorptionsCount,3);
    mRGCsNum = size(obj.rgcRFcenterConeConnectivityMatrix,2);
    
    % Allocate memory for the responses
    responses = zeros(nTrials, nTimePoints, mRGCsNum);

    % Compute the response of each mRGC
    parfor iRGC = 1:mRGCsNum
        connectivityVector = full(squeeze(obj.rgcRFcenterConeConnectivityMatrix(:, iRGC)));
        inputConeIndices = find(connectivityVector > 0.0001);
        inputConeWeights = reshape(connectivityVector(inputConeIndices), [1 1 numel(inputConeIndices)]);
        % Sum the weighted cone responses
        tmp = sum(bsxfun(@times, noiseFreeAbsorptionsCount(:,:,inputConeIndices), inputConeWeights),3);
        responses(:,:,iRGC) = tmp;
    end

end
