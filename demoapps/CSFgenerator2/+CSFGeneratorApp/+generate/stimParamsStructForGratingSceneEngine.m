function [sParams, coneMosaicIsTooSmall] = stimParamsStructForGratingSceneEngine(app, currentSpatialFrequency)
    % Form sParams struct
    sParams = struct(...
                'contrast', 1.0, ...
                'chromaDir', [app.stimParams.LconeContrast app.stimParams.MconeContrast app.stimParams.SconeContrast]/100, ...
                'meanLuminanceCdPerM2', app.stimParams.meanLuminanceCdM2, ...
                'sf', currentSpatialFrequency, ...
                'fovDegs', app.stimParams.sizeDegs, ...
                'orientation', app.stimParams.orientationDegs, ...
                'spatialEnvelope', app.stimParams.spatialEnvelope, ...
                'spatialPhase', app.stimParams.spatialPhaseDegs, ...
                'pixelsNum', app.stimParams.resolutionPixels, ...
                'minPixelsNumPerCycle', app.stimParams.minPixelsNumPerCycle, ...
                'spectralSupport', (app.stimParams.wavelengthSupportMin : app.stimParams.wavelengthSupportStepSize: app.stimParams.wavelengthSupportMax), ...
                'presentationMode', app.stimParams.presentationMode, ...
                'duration', app.stimParams.durationSec, ...
                'warningInsteadOfErrorOnOutOfGamut', true ...
     );

     switch (app.stimParams.spatialEnvelope)
            case 'disk'
                sParams.spatialEnvelopeRadiusDegs = app.stimParams.sizeDegs/2;
            case 'rect'
                sParams.spatialEnvelopeRadiusDegs = app.stimParams.sizeDegs/2;
            case 'soft'
                % In this case, the stimulus is a Gabor, and spatialEnvelopeRadiusDegs represents 1 sigma.
                if (strcmp(app.csfParams.constantParameter, 'constant cycles'))
                    sigmaDegs = app.csfParams.numberOfConstantCycles*1/currentSpatialFrequency;
                    bandwidthOctaves = bandwidthOctavesFromSigma(sigmaDegs , currentSpatialFrequency);
                elseif (strcmp(app.csfParams.constantParameter, 'constant size'))
                    sigmaDegs = app.stimParams.sizeDegs/6;
                end
                sParams.spatialEnvelopeRadiusDegs = sigmaDegs; 
            otherwise
                error('Unknown spatial envelope: ''%s''.', app.stimParams.spatialEnvelope);
     end
     
     % Check whether the cone mosaic is large enough for the stimulus
     coneMosaicIsTooSmall = false;
     maxStimSize = CSFGeneratorApp.compute.stimulusSizeAtLowestSpatialFrequency(app);
     if (min(app.coneMosaicParams.sizeDegs) < maxStimSize)
        coneMosaicIsTooSmall = true;
     end
            
end

function sigma = sigmaFromBandwidthOctaves(b, sf)
    sigma = 1.0/(sf *pi) * sqrt(log(2)/2) * (2^b+1)/(2^b-1);
end

function b = bandwidthOctavesFromSigma(s, sf)
    n = s*sf*pi+sqrt(log(2)/2);
    d = s*sf*pi-sqrt(log(2)/2);
    b = log2(n/d);
end

