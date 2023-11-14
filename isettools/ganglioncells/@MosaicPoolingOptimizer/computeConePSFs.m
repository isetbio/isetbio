function [theLconePSF, theMconePSF, theSconePSF, spatialSupportDegs] = ...
    computeConePSFs(theRGCMosaic, opticsToEmploy, stimSizeDegs, stimXYpositionDegs, varargin)

    p = inputParser;
    p.addParameter('coneFundamentalsOptimizedForStimPosition', true, @islogical);
    p.parse(varargin{:});
    coneFundamentalsOptimizedForStimPosition = p.Results.coneFundamentalsOptimizedForStimPosition;

    % Retrieve the input cone mosaic
    theConeMosaic = theRGCMosaic.inputConeMosaic;

    % Retrieve the optics to employ
    switch (opticsToEmploy)
        case {'native', 'adaptive optics'}
             % Retrieve the native optics
             theOI = theRGCMosaic.theNativeOptics;
        case 'custom'
            % Retrieve the custom optics
             theOI = theRGCMosaic.theCustomOptics;

        otherwise
             error('Unknown optics: ''%s''.', opticsToEmploy);
    end


    if (coneFundamentalsOptimizedForStimPosition)
        % Compute custom cone fundamentals
        maxConesNumForAveraging = 3;
        customConeFundamentals = MosaicPoolingOptimizer.coneFundamentalsAtTargetPositionWithinConeMosaic(...
            theConeMosaic, theOI, stimXYpositionDegs, stimSizeDegs, maxConesNumForAveraging);
        wavelengthSupport = customConeFundamentals.wavelengthSupport;
        coneFundamentals = customConeFundamentals.quantalExcitationSpectra;
    
    else
        wavelengthSupport = theConeMosaic.wave;
        coneFundamentals = ieReadSpectra(fullfile(isetbioDataPath,'human','stockman'), displayWavelengths);
    end

    theLconePSF = [];
    theMconePSF = [];
    theSconePSF = [];
    theOptics = oiGet(theOI, 'optics');

    psfSupportMicrons = opticsGet(theOptics,'psf support','um');
    micronsPerDegree = oiGet(theOI,'optics dist per deg')*10^6; % um/deg
    xGridDegs = psfSupportMicrons{1}/micronsPerDegree;
    yGridDegs = psfSupportMicrons{2}/micronsPerDegree;
    spatialSupportDegs = xGridDegs(1,:);
    
    for iWave = 1:numel(wavelengthSupport)
        targetWavelength = wavelengthSupport(iWave);
        psfAtTargetWavelength = opticsGet(theOptics,'psf data',targetWavelength);

        if (iWave == 1)     
             theLconePSF = psfAtTargetWavelength * coneFundamentals(iWave,1);
             theMconePSF = psfAtTargetWavelength * coneFundamentals(iWave,2);
             theSconePSF = psfAtTargetWavelength * coneFundamentals(iWave,3);
        else
            theLconePSF = theLconePSF + psfAtTargetWavelength * coneFundamentals(iWave,1);
            theMconePSF = theMconePSF + psfAtTargetWavelength * coneFundamentals(iWave,2);
            theSconePSF = theSconePSF + psfAtTargetWavelength * coneFundamentals(iWave,3);
        end
    end

    theLconePSF = flipud(theLconePSF);
    theMconePSF = flipud(theMconePSF);
    theSconePSF = flipud(theSconePSF);

end
