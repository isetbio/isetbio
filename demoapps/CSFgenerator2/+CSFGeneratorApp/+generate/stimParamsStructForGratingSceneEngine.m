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
 
     spatialPeriodDegs = 1.0/currentSpatialFrequency;
     if (strcmp(app.csfParams.constantParameter, 'constant cycles'))
        stimSizeDegs = app.csfParams.numberOfConstantCycles * spatialPeriodDegs;
     else
        stimSizeDegs = app.stimParams.sizeDegs;
     end
        
     switch (app.stimParams.spatialEnvelope)
            case 'disk'
                sParams.spatialEnvelopeRadiusDegs = stimSizeDegs/2;
            case 'rect'
                sParams.spatialEnvelopeRadiusDegs = stimSizeDegs/2;
            case 'soft'
                % In this case, where the stimulus is a Gabor,
                % spatialEnvelopeRadiusDegs represents 1 sigma.
                
                if (strcmp(app.csfParams.constantParameter, 'constant cycles'))
                    bandwidthOctaves = 2/app.csfParams.numberOfConstantCycles;
                else
                    bandwidthOctaves = 2/(stimSizeDegs * currentSpatialFrequency);
                end
                sigma = sqrt(log(2))*(2^bandwidthOctaves + 1)/(sqrt(2)*pi*currentSpatialFrequency*(2^bandwidthOctaves-1));
                sigmaCycles = sigma / spatialPeriodDegs;
                sParams.spatialEnvelopeRadiusDegs = sigma; 
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

