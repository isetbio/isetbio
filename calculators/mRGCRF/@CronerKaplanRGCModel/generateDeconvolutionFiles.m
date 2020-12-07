function generateDeconvolutionFiles(obj, deconvolutionOpticsParams, subregion, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('examinedConesNumInRFCenter', 1:30);
    p.addParameter('sensitivityRangeOverWhichToMatchSFtuning', [0.9 0.05]);
    p.addParameter('visualizeFits', false);
    p.addParameter('exportFig', false);
    p.parse(varargin{:});

    % Validate the deconvolutionOpticsParams
    obj.validateDeconvolutionOpticsParams(deconvolutionOpticsParams);
    
    examinedConesNumInRFCenter = p.Results.examinedConesNumInRFCenter;
    sensitivityRangeOverWhichToMatchSFtuning = p.Results.sensitivityRangeOverWhichToMatchSFtuning;
    visualizeFits = p.Results.visualizeFits;
    exportFig = p.Results.exportFig;
    imposedRefractionErrorDiopters = 0;
    pupilDiamMM = 3.0;
    
    deconvolutionEccentricities = obj.deconvolutionEccs;
    
    if (strcmp(subregion, 'center'))
        % Parfor on the eccentricities
        parfor eccIndex = 1:numel(deconvolutionEccentricities)
            patchEccRadiusDegs = -deconvolutionEccentricities(eccIndex);
            generateDeconvolutionFilesForTheCenter(obj, patchEccRadiusDegs, examinedConesNumInRFCenter, ...
                sensitivityRangeOverWhichToMatchSFtuning, deconvolutionOpticsParams, imposedRefractionErrorDiopters, pupilDiamMM, visualizeFits, exportFig);
        end
    else
        % No parfor on the eccentricities. Parfor is called within generateDeconvolutionFilesForTheSurround
        for eccIndex = 1:numel(deconvolutionEccentricities)
            patchEccRadiusDegs = -deconvolutionEccentricities(eccIndex);
            generateDeconvolutionFilesForTheSurround(obj, patchEccRadiusDegs, examinedConesNumInRFCenter, ...
                deconvolutionOpticsParams, imposedRefractionErrorDiopters, pupilDiamMM, visualizeFits, exportFig);
        end
    end
    
end

function generateDeconvolutionFilesForTheCenter(obj, patchEccRadiusDegs, examinedConesNumInRFCenter, ...
    sensitivityRangeOverWhichToMatchSFtuning, deconvolutionOpticsParams, imposedRefractionErrorDiopters, pupilDiamMM, visualizeFits, exportFig)
    
    % Extra deconvolution optics params
    subjectIDs = deconvolutionOpticsParams.PolansWavefrontAberrationSubjectIDsToCompute;
    quadrants = deconvolutionOpticsParams.quadrantsToCompute;
    subjectsNum = numel(subjectIDs);
    quadrantsNum = numel(quadrants);
   
    resetPlotLabOnExit = false;
    
    deconvolutionStruct = cell(quadrantsNum, subjectsNum);
    
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
        extraPatchSizeDegs = 0;
        [conePosDegs, coneAperturesDegs, micronsPerDegree, wavelengthSampling] = ...
            generateConePositionsForPatchAtEccentricity(patchEccDegs, extraPatchSizeDegs, 'coneMosaicResamplingFactor', 1);
    
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
            deconvolutionStruct{qIndex, sIndex} = ...
                obj.performDeconvolutionAnalysisForRFcenter(examinedConesNumInRFCenter, ...
                sensitivityRangeOverWhichToMatchSFtuning, conePosDegs, coneAperturesDegs, thePSF, thePSFsupportDegs, visualizeFits, exportFig, ...
                quadrants{qIndex}, PolansSubjectID, patchEccRadiusDegs);
        end % sIndex
    end % qIndex
    
    % Save decolvolution struct for the center
    dataFileName = obj.deconvolutionDataFileName(patchEccRadiusDegs, imposedRefractionErrorDiopters, 'center');
    fprintf('Saving data to %s\n', dataFileName);
    save(dataFileName, ...
            'deconvolutionStruct', ...
            'subjectIDs', 'quadrants');  
        
    if (resetPlotLabOnExit)
        obj.setupPlotLab(-1);
    end
    
end


function generateDeconvolutionFilesForTheSurround(obj, patchEccRadiusDegs, examinedConesNumInRFCenter, ...
    deconvolutionOpticsParams, imposedRefractionErrorDiopters, pupilDiamMM, visualizeFits, exportFig)

    % Extra deconvolution optics params
    subjectIDs = deconvolutionOpticsParams.PolansWavefrontAberrationSubjectIDsToCompute;
    quadrants = deconvolutionOpticsParams.quadrantsToCompute;
    subjectsNum = numel(subjectIDs);
    quadrantsNum = numel(quadrants);
    
    % Load deconvolutionStruct for the center
    dataFileName = obj.deconvolutionDataFileName(patchEccRadiusDegs, imposedRefractionErrorDiopters, 'center');
    fprintf('Loading decolvolution data from the center %s\n', dataFileName);
    load(dataFileName, 'deconvolutionStruct');  
    deconvolutionStructForTheCenter = deconvolutionStruct;
    clear 'deconvolutionStruct'
    
    deconvolutionStruct = cell(quadrantsNum, subjectsNum);
    
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
        extraPatchSizeDegs = 0.5*sqrt(max(examinedConesNumInRFCenter)) * 2 * 3 * obj.surroundRadiusFunction(obj.surroundRadiusParams, sqrt(sum(patchEccDegs.^2,2)));
        
        [conePosDegs, coneAperturesDegs, micronsPerDegree, wavelengthSampling] = ...
            generateConePositionsForPatchAtEccentricity(patchEccDegs, extraPatchSizeDegs, 'coneMosaicResamplingFactor', 1);
    
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
                deconvolutionStructForTheCenter{qIndex, sIndex}, examinedConesNumInRFCenter, ...
                conePosDegs, coneAperturesDegs, thePSF, thePSFsupportDegs, visualizeFits, exportFig, ...
                quadrants{qIndex}, PolansSubjectID, patchEccRadiusDegs);
            
        end % sIndex
    end % qIndex
    
    % Save deconvolution struct for the surround
    dataFileName = obj.deconvolutionDataFileName(patchEccRadiusDegs, imposedRefractionErrorDiopters, 'surround');
    fprintf('Saving data to %s\n', dataFileName);
    save(dataFileName, ...
            'deconvolutionStruct', ...
            'subjectIDs', 'quadrants');  
        
    if (resetPlotLabOnExit)
        obj.setupPlotLab(-1);
    end
    
end


function [conePosDegs, coneAperturesDegs, micronsPerDegree, wavelengthSampling] = generateConePositionsForPatchAtEccentricity(...
    patchEccDegs, extraPatchSizeDegs, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('coneMosaicResamplingFactor', 9, @isnumeric);
    p.parse(varargin{:});
    
    patchSizeDegs = extraPatchSizeDegs*[1 1]; 
    [theConeMosaic, theConeMosaicMetaData] = CronerKaplanRGCModel.generateConeMosaicForDeconvolution(...
                patchEccDegs, patchSizeDegs, ...
                'coneMosaicResamplingFactor', p.Results.coneMosaicResamplingFactor, ...
                'sizeUnits', 'degrees', ...
                'mosaicGeometry', 'regular');
            
    % Subtract the patchCenter because the PSF is also going to be centered at (0,0) for this analysis
    conePosDegs = bsxfun(@minus, theConeMosaicMetaData.conePositionsDegs, theConeMosaicMetaData.coneMosaicEccDegs);
    % Make sure central-most cone is at (0,0)
    [~,idx] = min(sqrt(sum(conePosDegs.^2,2)));
    conePosDegs = bsxfun(@minus, conePosDegs, conePosDegs(idx,:));
    coneAperturesDegs = theConeMosaicMetaData.coneAperturesDegs;

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