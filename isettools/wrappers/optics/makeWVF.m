function theWVF = makeWVF(wavefrontSpatialSamples, zcoeffsMicrons, measWavelength, wavelengthsToCompute, ...
    measPupilDiameterMM, calcPupilDiameterMM, umPerDegree, name, varargin)

% The only place this is called from is computePSFandOTF
%
    % Parse input
    p = inputParser;
    p.addParameter('flipPSFUpsideDown', false, @islogical);
    p.addParameter('rotatePSF90degs', false, @islogical);
    p.addParameter('upsampleFactor', [], @(x)(isempty(x) || ((isnumeric(x))&&(numel(x)==1)&&(x>0))));
    p.addParameter('humanlca',true,@islogical);
    p.parse(varargin{:});
    flipPSFUpsideDown = p.Results.flipPSFUpsideDown;
    rotatePSF90degs = p.Results.rotatePSF90degs;
    upsampleFactor = p.Results.upsampleFactor;
    LCA = p.Results.humanlca;


    % Hack to change the psfsamplespacing according to the upsample factor
    referencePupilPlaneSize = 16.212;
    if (~isempty(upsampleFactor))
        wavefrontSpatialSamples = wavefrontSpatialSamples*upsampleFactor;
        referencePupilPlaneSize = referencePupilPlaneSize * upsampleFactor;
    end


    theWVF = wvfCreate(...
    			'umPerDegree', umPerDegree, ...
                'calc wavelengths',wavelengthsToCompute,...
                'measuredpupil', measPupilDiameterMM, ...
                'calc pupil size',calcPupilDiameterMM, ...
                'spatialsamples', wavefrontSpatialSamples, ...
                'ref pupil plane size', referencePupilPlaneSize, ...
                'zcoeffs', zcoeffsMicrons,...
                'measured wl', measWavelength, ...
                'name', name, ...
                'flipPSFUpsideDown', flipPSFUpsideDown, ...
                'rotatePSF90degs', rotatePSF90degs);

    if (~isempty(upsampleFactor))
        arcminPerSample = wvfGet(theWVF,'psf angle per sample','min',measWavelength);
        fprintf('PSF angular resolution: %f min\n', arcminPerSample);
    end
    
    % Now compute the PSF.  Also set customLca field if human LCA is on.
    if (LCA)
        theWVF = wvfSet(theWVF,'lcaMethod','human');
    else
        theWVF = wvfSet(theWVF,'lcaMethod','none');
    end
    theWVF = wvfCompute(theWVF);
end