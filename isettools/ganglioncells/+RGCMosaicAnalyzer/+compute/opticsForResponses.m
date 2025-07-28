function [theOI, thePSF] = opticsForResponses(theMRGCMosaic, whichOptics, customRefractionDiopters, visualizsPSFonTopOfConeMosaic)

	switch (whichOptics)
        case 'nativeOptics'
            % Generate the optics that where used to optimize the mosaic
            [theOI, thePSF] = theMRGCMosaic.nativeOI(...
                'visualizePSF', true);

        case 'adaptiveOptics6MM'
            [theOI, thePSF] = theMRGCMosaic.nativeOI(...
                'opticsModification', whichOptics, ...
                'visualizePSF', true);

        case 'adaptiveOptics6MMwithLCA'
            [theOI, thePSF] = theMRGCMosaic.nativeOI(...
                'opticsModification', whichOptics, ...
                'visualizePSF', true);

       	case 'customRefraction'
            [theOI, thePSF] = theMRGCMosaic.nativeOI(...
                'opticsModification', whichOptics, ...
                'customRefractionDiopters', customRefractionDiopters, ...
                'visualizePSF', true, ...
                'visualizedWavelengths', 450:20:650);

        case 'refractionResidualWithRespectToNativeOptics'
            % Retrieve the defocus required to optimized the Strehl ratio
           [~,~,theOptimalStrehlRatioDefocusDiopters] = theMRGCMosaic.nativeOI(...
                'visualizePSF', true);
            theOptimalStrehlRatioDefocusDiopters
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
                'visualizePSF', true);

   	end % switch

   	if (visualizsPSFonTopOfConeMosaic)

        % Compute the 2-deg  Stockman cone fundamentals weighted PSFs
        [theLconeWeightedPSF, theMconeWeightedPSF, theSconeWeightedPSF] = ...
            RGCMosaicConstructor.helper.optics.coneWeightedPSFs(thePSF);

        RGCMosaicAnalyzer.visualize.coneWeightedPSFonTopOfConeMosaic(1, theMRGCMosaic.inputConeMosaic, theLconeWeightedPSF, sprintf('%s_LconeWeightedPSF.pdf',whichOptics), 'L-cone weighted PSF');
        RGCMosaicAnalyzer.visualize.coneWeightedPSFonTopOfConeMosaic(2, theMRGCMosaic.inputConeMosaic, theMconeWeightedPSF, sprintf('%s_MconeWeightedPSF.pdf',whichOptics), 'M-cone weighted PSF');
        RGCMosaicAnalyzer.visualize.coneWeightedPSFonTopOfConeMosaic(3, theMRGCMosaic.inputConeMosaic, theSconeWeightedPSF, sprintf('%s_SconeWeightedPSF.pdf',whichOptics), 'S-cone weighted PSF');
    end

end
