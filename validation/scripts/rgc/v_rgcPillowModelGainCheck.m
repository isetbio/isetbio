% v_rgcNumericalCheck
% 
% A validation for the precomputed retinal ganglion cell response models.
% 
% TRG 10/2017, (c) isetbio team

%% Initialize
clear; close all; ieInit;

%% Define cell types and sample rates to check
cellTypes = {'onmidget', 'onparasol', 'onsbc', 'offmidget', 'offparasol'};
sampleRates = [1/100, 1/1000];

%% Check the gain of the Pillow rgc temporal response model for different sampling rates
for i=1:length(cellTypes)
    for j=1:length(sampleRates)
        params = {};
        params.filterDuration = 0.4;
        params.samplingTime = sampleRates(j);
        params.cellType = cellTypes{i};
        [rgcFilter,rgcTime ] = rgcImpulseResponsePillow(params);
        rgcGain = sum(rgcFilter(:));

        if abs(rgcGain - 1.0) < 1e-6
            fprintf('Pillow %s filter with sample rate %0.4f has correct unit gain.\n', cellTypes{i}, sampleRates(j));
        else
            fprintf('Pillow %s filter with sample rate %0.4f has gain of %0.4f which should be 1.0\n', cellTypes{i}, sampleRates(j), rgcGain);
        end
    end
end