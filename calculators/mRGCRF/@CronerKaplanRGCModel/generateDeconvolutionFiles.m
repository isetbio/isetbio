function generateDeconvolutionFiles(obj, deconvolutionOpticsParams, varargin)

    defaultEccTested = [0 0.25 0.5 1 1.5 2.0 2.5 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25];
    defaultCenterConesNumTested = 1:100;
    
    % Parse input
    p = inputParser;
    p.addParameter('eccTested', defaultEccTested);
    p.addParameter('conesNumInRFcenterTested', defaultCenterConesNumTested);
            
    p.parse(varargin{:});

    % Validate the deconvolutionOpticsParams
    obj.validateDeconvolutionOpticsParams(deconvolutionOpticsParams);
    
    eccTested = -(p.Results.eccTested);
    conesNumInRFcenterTested = p.Results.conesNumInRFcenterTested;
    
    imposedRefractionErrorDiopters = 0;
    pupilDiamMM = 3.0;
    
    for eccIndex = 1:numel(eccTested)
        ecc = eccTested(eccIndex);
        doIt(obj, ecc, conesNumInRFcenterTested , deconvolutionOpticsParams, imposedRefractionErrorDiopters, pupilDiamMM);
    end
    
end

function doIt(obj, patchEccDegs, conesNumInRFcenterTested, deconvolutionOpticsParams, imposedRefractionErrorDiopters, pupilDiamMM)
    
    % Extra deconvolution optics params
    subjectIDs = deconvolutionOpticsParams.PolansWavefrontAberrationSubjectIDsToCompute;
    quadrants = deconvolutionOpticsParams.quadrantsToCompute;
    subjectsNum = numel(subjectIDs);
    quadrantsNum = numel(quadrants);
   
    summaryFigureNo = round(1000 + abs(patchEccDegs));
    hSummaryFig = figure(summaryFigureNo); clf;
    resetPlotLabOnExit = false;
    
    for qIndex = 1:quadrantsNum

        switch (quadrants{qIndex})
            case 'horizontal'
                % horizontal meridian (nasal)
                eccXrange = patchEccDegs(1)*[1 1];
                eccYrange = 0*[1 1];
            case 'superior'
                % vertical meridian (superior)
                eccYrange = patchEccDegs(1)*[1 1];
                eccXrange = 0*[1 1];
            case 'inferior'
                % vertical meridian (inferior)
                eccYrange = -patchEccDegs(1)*[1 1];
                eccXrange = 0*[1 1];
            otherwise
                error('Unknown Polans quadrant: ''%s''.', eccQuadrant);
        end
        
        % Generate cone positions appropriate for the eccentricity
        patchEccDegs = [eccXrange(1), eccYrange(1)];
        [conePosDegs, coneAperturesDegs, micronsPerDegree, wavelengthSampling] = generateConePositionsForPatchAtEccentricity(...
            patchEccDegs, 'coneMosaicResamplingFactor', 1);
    
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
                CronerKaplanRGCModel.performDeconvolutionAnalysis(conesNumInRFcenterTested, ...
                conePosDegs, coneAperturesDegs, thePSF, thePSFsupportDegs);
                  
                
%                 if (visualizeAnalysis)
%                     % Extract profile at peak
%                     row = round(size(rfPoolingInRetinalSpace,1)/2);
%                     rfPoolingRetinalSpaceProfile = squeeze(rfPoolingInRetinalSpace(row,:));
%                     rfPoolingVisualSpaceProfile = squeeze(rfPoolingInVisualSpace(row,:));
%                     maxX = 0.1;
%                     visualizeFit(thePSFsupportDegs, thePSF, subjectID, eccXrange, eccYrange, ...
%                         rfPoolingInRetinalSpaceNorm, rfPoolingInVisualSpaceNorm, ...
%                         rfPoolingRetinalSpaceProfile, rfPoolingVisualSpaceProfile, ...
%                         retinalRadius(retinalRadiusIndex,qIndex, sIndex), visualRadius(retinalRadiusIndex,qIndex, sIndex), ...
%                         visualGain(retinalRadiusIndex,qIndex, sIndex), ellipseInRetinalSpace, ellipseInVisualSpace, cMap, maxX);
%                 end
        end % sIndex
    end % qIndex
    
    % Save data
    dataFileName = fullfile(obj.psfDeconvolutionDir, sprintf('ecc_%2.1f_%2.1f_centerDeconvolutions_refractionError_%2.2fD.mat', eccXrange(1), eccYrange(1), imposedRefractionErrorDiopters));
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

