function [theOI, thePSF] = opticsForResponses(theMRGCMosaic, whichOptics, customRefractionDiopters, visualizePSFonTopOfConeMosaic, varargin)

    p = inputParser;
    % Optional params
    p.addParameter('wavefrontSpatialSamples', [], @(x)(isempty(x)||isnumeric(x)));
    p.parse(varargin{:});
    wavefrontSpatialSamples = p.Results.wavefrontSpatialSamples;

	switch (whichOptics)

        case 'nativeOpticsNoStrehlRatioOptimization'
            % Generate the optics that where used to optimize the mosaic
            [theOI, thePSF] = theMRGCMosaic.nativeOI(...
                'opticsModification', whichOptics, ...
                'wavefrontSpatialSamples', wavefrontSpatialSamples, ...
                'visualizePSF', visualizePSFonTopOfConeMosaic);

        case 'nativeOptics'
            % Generate the optics that where used to optimize the mosaic
            fprintf('\n\nComputing the defocus value that maximizes the Strehl ratio. Please wait ...\n\n')
            [theOI, thePSF, theOptimalStrehlRatioDefocusDiopters] = theMRGCMosaic.nativeOI(...
                'wavefrontSpatialSamples', wavefrontSpatialSamples, ...
                'visualizePSF', visualizePSFonTopOfConeMosaic);

            if (isempty(theOptimalStrehlRatioDefocusDiopters))
                error('The native optics in the compute-ready mosaic were not constructed so as to optimize the Strehl ratio.')
            end

        case 'adaptiveOptics6MM'
            [theOI, thePSF] = theMRGCMosaic.nativeOI(...
                'opticsModification', whichOptics, ...
                'wavefrontSpatialSamples', wavefrontSpatialSamples, ...
                'visualizePSF', visualizePSFonTopOfConeMosaic);

        case 'adaptiveOptics6MMwithLCA'
            [theOI, thePSF] = theMRGCMosaic.nativeOI(...
                'opticsModification', whichOptics, ...
                'wavefrontSpatialSamples', wavefrontSpatialSamples, ...
                'visualizePSF', visualizePSFonTopOfConeMosaic);

       	case 'customRefraction'
            [theOI, thePSF] = theMRGCMosaic.nativeOI(...
                'opticsModification', whichOptics, ...
                'customRefractionDiopters', customRefractionDiopters, ...
                'wavefrontSpatialSamples', wavefrontSpatialSamples, ...
                'visualizePSF', visualizePSFonTopOfConeMosaic, ...
                'visualizedWavelengths', 450:20:650);

        case {'refractionResidualWithRespectToNativeOptics', 'loadComputeReadyRGCMosaic'}
            % Retrieve the defocus required to optimized the Strehl ratio
            fprintf('\n\nComputing the defocus value that maximizes the Strehl ratio. Please wait ...\n\n')
            [~,~,theOptimalStrehlRatioDefocusDiopters] = theMRGCMosaic.nativeOI(...
                'wavefrontSpatialSamples', wavefrontSpatialSamples, ...
                'visualizePSF', visualizePSFonTopOfConeMosaic);

            if (isempty(theOptimalStrehlRatioDefocusDiopters))
                error('The native optics in the compute-ready mosaic were not constructed so as to optimize the Strehl ratio.')
            end

           % Add in custom refraction (relative to theOptimalStrehlRatioDefocusDiopters)
           customRefractionDiopters = ...
                theOptimalStrehlRatioDefocusDiopters + ...
                customRefractionDiopters;

           [theOI, thePSF] = theMRGCMosaic.nativeOI(...
                'opticsModification', 'customRefraction', ...
                'customRefractionDiopters', customRefractionDiopters, ...
                'wavefrontSpatialSamples', wavefrontSpatialSamples, ...
                'visualizePSF', visualizePSFonTopOfConeMosaic);

        otherwise
            error('Not allowed to generate optics for responses with this specification: ''%s''.', whichOptics);
   	end % switch

   	if (visualizePSFonTopOfConeMosaic)

        % Compute the 2-deg  Stockman cone fundamentals weighted PSFs
        [theLconeWeightedPSF, theMconeWeightedPSF, theSconeWeightedPSF] = ...
            RGCMosaicConstructor.helper.optics.coneWeightedPSFs(thePSF);

        RGCMosaicAnalyzer.visualize.coneWeightedPSFonTopOfConeMosaic(1, theMRGCMosaic.inputConeMosaic, theLconeWeightedPSF, sprintf('%s_LconeWeightedPSF.pdf',whichOptics), 'L-cone weighted PSF');
        RGCMosaicAnalyzer.visualize.coneWeightedPSFonTopOfConeMosaic(2, theMRGCMosaic.inputConeMosaic, theMconeWeightedPSF, sprintf('%s_MconeWeightedPSF.pdf',whichOptics), 'M-cone weighted PSF');
        RGCMosaicAnalyzer.visualize.coneWeightedPSFonTopOfConeMosaic(3, theMRGCMosaic.inputConeMosaic, theSconeWeightedPSF, sprintf('%s_SconeWeightedPSF.pdf',whichOptics), 'S-cone weighted PSF');
    end

end
