function generateDeconvolutionFiles(obj, deconvolutionOpticsParams, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('conesNumInRFcenterTested', 1:30);
    p.addParameter('visualizeFits', false);
    p.parse(varargin{:});

    % Validate the deconvolutionOpticsParams
    obj.validateDeconvolutionOpticsParams(deconvolutionOpticsParams);
    
    deconvolutionEccs = -obj.deconvolutionEccs;
    conesNumInRFcenterTested = p.Results.conesNumInRFcenterTested;
    visualizeFits = p.Results.visualizeFits;
    
    imposedRefractionErrorDiopters = 0;
    pupilDiamMM = 3.0;
    
    for eccIndex = 1:numel(deconvolutionEccs)
        patchEccRadiusDegs = deconvolutionEccs(eccIndex);
        generateDeconvolutionFilesForTheCenter(obj, patchEccRadiusDegs, conesNumInRFcenterTested, deconvolutionOpticsParams, imposedRefractionErrorDiopters, pupilDiamMM, visualizeFits);
        generateDeconvolutionFilesForTheSurround(obj, patchEccRadiusDegs, deconvolutionOpticsParams, imposedRefractionErrorDiopters, pupilDiamMM, visualizeFits);
    end
    
end

function generateDeconvolutionFilesForTheCenter(obj, patchEccRadiusDegs, conesNumInRFcenterTested, deconvolutionOpticsParams, imposedRefractionErrorDiopters, pupilDiamMM, visualizeFits)
    
    % Extra deconvolution optics params
    subjectIDs = deconvolutionOpticsParams.PolansWavefrontAberrationSubjectIDsToCompute;
    quadrants = deconvolutionOpticsParams.quadrantsToCompute;
    subjectsNum = numel(subjectIDs);
    quadrantsNum = numel(quadrants);
   
    resetPlotLabOnExit = false;
    
    for qIndex = 1:quadrantsNum

        switch (quadrants{qIndex})
            case 'horizontal'
                % horizontal meridian (nasal)
                eccXrange = patchEccRadiusDegs*[1 1];
                eccYrange = 0*[1 1];
            case 'superior'
                % vertical meridian (superior)
                eccYrange = patchEccRadiusDegs*[1 1];
                eccXrange = 0*[1 1];
            case 'inferior'
                % vertical meridian (inferior)
                eccYrange = -patchEccRadiusDegs*[1 1];
                eccXrange = 0*[1 1];
            otherwise
                error('Unknown Polans quadrant: ''%s''.', eccQuadrant);
        end
        
        % Generate cone positions appropriate for the eccentricity
        patchEccDegs = [eccXrange(1), eccYrange(1)];
        [conePosDegs, coneAperturesDegs, micronsPerDegree, wavelengthSampling] = ...
            generateConePositionsForPatchAtEccentricity(patchEccDegs, 'coneMosaicResamplingFactor', 1);
    
        for sIndex = 1:subjectsNum
            % Compute the Polans subject PSF at the desired eccentricity
            PolansSubjectID = deconvolutionOpticsParams.PolansWavefrontAberrationSubjectIDsToCompute(sIndex);
            
            [thePSF, thePSFsupportDegs] = generatePSFForPatchAtEccentricity(PolansSubjectID, ...
                imposedRefractionErrorDiopters, pupilDiamMM, wavelengthSampling, ...
                micronsPerDegree, patchEccDegs, ...
                'centerPSFatZero', true);
            
            if (qIndex == 1) && (sIndex == 1)
                obj.setupPlotLab(0, 20, 12);
                resetPlotLabOnExit = true;
            end
            
            % Convolve different retinal pooling regions and compute the visually-mapped pooling region
            deconvolutionStruct{ qIndex, sIndex} = ...
                obj.performDeconvolutionAnalysisForRFcenter(conesNumInRFcenterTested, ...
                conePosDegs, coneAperturesDegs, thePSF, thePSFsupportDegs, visualizeFits);
        end % sIndex
    end % qIndex
    
    % Save decolvolution struct for the center
    dataFileName = fullfile(obj.psfDeconvolutionDir, sprintf('ecc_%2.1f_centerDeconvolutions_refractionError_%2.2fD.mat', patchEccRadiusDegs, imposedRefractionErrorDiopters));
    fprintf('Saving data to %s\n', dataFileName);
    save(dataFileName, ...
            'deconvolutionStruct', ...
            'subjectIDs', 'quadrants');  
        
    if (resetPlotLabOnExit)
        obj.setupPlotLab(-1);
    end
    
end


function generateDeconvolutionFilesForTheSurround(obj, patchEccRadiusDegs, deconvolutionOpticsParams, imposedRefractionErrorDiopters, pupilDiamMM, visualizeFits)
    % Extra deconvolution optics params
    subjectIDs = deconvolutionOpticsParams.PolansWavefrontAberrationSubjectIDsToCompute;
    quadrants = deconvolutionOpticsParams.quadrantsToCompute;
    subjectsNum = numel(subjectIDs);
    quadrantsNum = numel(quadrants);
    
    % Load deconvolutionStruct for the center
    dataFileName = fullfile(obj.psfDeconvolutionDir, sprintf('ecc_%2.1f_centerDeconvolutions_refractionError_%2.2fD.mat', patchEccRadiusDegs, imposedRefractionErrorDiopters));
    fprintf('Loading decolvolution data from the center %s\n', dataFileName);
    load(dataFileName, 'deconvolutionStruct');  
    deconvolutionStructForTheCenter = deconvolutionStruct;
    
    for qIndex = 1:quadrantsNum

        switch (quadrants{qIndex})
            case 'horizontal'
                % horizontal meridian (nasal)
                eccXrange = patchEccRadiusDegs*[1 1];
                eccYrange = 0*[1 1];
            case 'superior'
                % vertical meridian (superior)
                eccYrange = patchEccRadiusDegs*[1 1];
                eccXrange = 0*[1 1];
            case 'inferior'
                % vertical meridian (inferior)
                eccYrange = -patchEccRadiusDegs*[1 1];
                eccXrange = 0*[1 1];
            otherwise
                error('Unknown Polans quadrant: ''%s''.', eccQuadrant);
        end
        
        % Generate cone positions appropriate for the eccentricity
        patchEccDegs = [eccXrange(1), eccYrange(1)];
        [conePosDegs, coneAperturesDegs, micronsPerDegree, wavelengthSampling] = ...
            generateConePositionsForPatchAtEccentricity(patchEccDegs, 'coneMosaicResamplingFactor', 1);
    
        for sIndex = 1:subjectsNum
            % Compute the Polans subject PSF at the desired eccentricity
            PolansSubjectID = deconvolutionOpticsParams.PolansWavefrontAberrationSubjectIDsToCompute(sIndex);
            
            [thePSF, thePSFsupportDegs] = generatePSFForPatchAtEccentricity(PolansSubjectID, ...
                imposedRefractionErrorDiopters, pupilDiamMM, wavelengthSampling, ...
                micronsPerDegree, patchEccDegs, ...
                'centerPSFatZero', true);
            
            if (qIndex == 1) && (sIndex == 1)
                obj.setupPlotLab(0, 20, 12);
                resetPlotLabOnExit = true;
            end
            
            % Convolve different retinal pooling regions and compute the visually-mapped pooling region
            deconvolutionStruct{qIndex, sIndex} = obj.performDeconvolutionAnalysisForRFsurround(...
                deconvolutionStructForTheCenter{qIndex, sIndex}, ...
                conePosDegs, coneAperturesDegs, thePSF, thePSFsupportDegs, visualizeFits);
            
        end % sIndex
    end % qIndex
    
    % Save deconvolution struct for the surround
    dataFileName = fullfile(obj.psfDeconvolutionDir, sprintf('ecc_%2.1f_surroundDeconvolutions_refractionError_%2.2fD.mat', patchEccRadiusDegs, imposedRefractionErrorDiopters));
    fprintf('Saving data to %s\n', dataFileName);
    save(dataFileName, ...
            'deconvolutionStruct', ...
            'subjectIDs', 'quadrants');  
        
    if (resetPlotLabOnExit)
        obj.setupPlotLab(-1);
    end
    
end


function [conePosDegs, coneAperturesDegs, micronsPerDegree, wavelengthSampling] = generateConePositionsForPatchAtEccentricity(patchEccDegs, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('coneMosaicResamplingFactor', 9, @isnumeric);
    p.parse(varargin{:});
    
    patchSizeDegs = [0 0];  % Just the max surround region as computed in mosaicsAndOpticsForEccentricity
    [theConeMosaic, theConeMosaicMetaData] = CronerKaplanRGCModel.generateConeMosaicForDeconvolution(...
                patchEccDegs, patchSizeDegs, ...
                'coneMosaicResamplingFactor', p.Results.coneMosaicResamplingFactor, ...
                'sizeUnits', 'degrees', ...
                'mosaicGeometry', 'regular');
            
    conePosDegs = WatsonRGCModel.rhoMMsToDegs(theConeMosaicMetaData.conePositionsMicrons*1e-3);
    % Subtract the patchCenter because the PSF is also going to be centered at (0,0) for this analysis
    conePosDegs = bsxfun(@minus, conePosDegs, theConeMosaicMetaData.coneMosaicEccDegs);
    % Make sure central-most cone is at (0,0)
    [~,idx] = min(sqrt(sum(conePosDegs.^2,2)));
    conePosDegs = bsxfun(@minus, conePosDegs, conePosDegs(idx,:));
    coneAperturesDegs = WatsonRGCModel.sizeRetinalMicronsToSizeDegs(theConeMosaicMetaData.coneAperturesMicrons, sqrt(sum((theConeMosaicMetaData.coneMosaicEccMicrons).^2,2)));

    micronsPerDegree = theConeMosaic.micronsPerDegree;
    wavelengthSampling = theConeMosaic.pigment.wave;
end

function [thePSFAtASingleWavelength, thePSFsupportDegs] = generatePSFForPatchAtEccentricity(PolansSubjectID, ...
    imposedRefractionErrorDiopters, pupilDiamMM, wavelengthSampling, micronsPerDegree, patchEccDegs, varargin)

     % Parse input
    p = inputParser;
    p.addParameter('centerPSFatZero', true, @islogical);
    p.parse(varargin{:});
    centerPSFatZero = p.Results.centerPSFatZero;
    
    [theOI, theFullPSF, thePSFsupportDegs] = CronerKaplanRGCModel.generatePolansOpticsForDeconcolution(...
        PolansSubjectID, imposedRefractionErrorDiopters, pupilDiamMM, ...
        wavelengthSampling, micronsPerDegree, patchEccDegs, ...
        'eccentricityUnits', 'degrees', ...
        'doNotZeroCenterPSF', ~centerPSFatZero);
    
    % Get PSF slice at target wavelength
    theOptics = oiGet(theOI, 'optics');
    targetWavelength = 550;
    thePSFAtASingleWavelength = opticsGet(theOptics,'psf data',targetWavelength);
end