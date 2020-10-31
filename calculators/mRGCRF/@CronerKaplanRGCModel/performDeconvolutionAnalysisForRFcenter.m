function deconvolutionStruct = performDeconvolutionAnalysisForRFcenter(obj, conesNumInRFcenterTested, ...
    sensitivityRangeOverWhichToMatchSFtuning, conePosDegs, coneAperturesDegs, thePSF, thePSFsupportDegs, visualizeFits, exportFig, ...
    quadrantName, subjectID, patchEccRadiusDegs)
    
    % Match the major spatial RF axis. This corresponding model will have larger RF centers, and therefore lower peak sensitivity
    spatialFrequencyAxisToMatch = 'major'; 
    
    % Match the minor spatial RF axis. This corresponding model will have smaller RF centers, and therefore higher peak sensitivity
    %spatialFrequencyAxisToMatch = 'minor';
    
    % Match the average of the minor & major spatial RF axis.
    %spatialFrequencyAxisToMatch = 'average';     

    fprintf('>>> RF center deconvolution: estimating visual characteristic radius by matching the **%s** spatial frequency axis.\n', spatialFrequencyAxisToMatch);
    
    % Flag indicating whether to overlay the derived matching Gaussian
    % profile on the retinal and visual cone images
    overlayMatchingGaussianProfileOnConeImages = true;
            
    % Interpolate PSF by a factor of 3
    upsampleFactor = 3;
    % Zero pad with a 0.0 degs margin on each size
    paddingMarginDegs = 0.0;
    [thePSFHR, thePSFsupportDegsHR] = ...
        CronerKaplanRGCModel.interpolatePSF(thePSF, thePSFsupportDegs, upsampleFactor, paddingMarginDegs);
    clear 'thePSFsupportDegs'
    
    % Generate the cone aperture profile
    [coneApertureProfile, ~] = obj.generateConeApertureProfileForDeconvolution(thePSFsupportDegsHR, coneAperturesDegs);
    
    % Examine a range of possible spatial arrangements of conesNumInRFcenterTested
    rfCenterPoolingSchemes = generateCenterPoolingSchemes(conesNumInRFcenterTested, coneAperturesDegs, conePosDegs);
    
    % DeconvolutionData is a container indexed by the number of cones in the RF center
    deconvolutionStruct.data = containers.Map();
    deconvolutionStruct.metaData = struct(...
        'patchEccRadiusDegs', patchEccRadiusDegs, ...
        'quadrantName',quadrantName, ...
        'subjectID', subjectID, ...
        'sensitivityRangeOverWhichToMatchSFtuning', sensitivityRangeOverWhichToMatchSFtuning ...
        );
    
    for poolingSchemeIndex = 1:numel(rfCenterPoolingSchemes)
        % Retrieve pooling scheme info
        conesNumInRFcenter = rfCenterPoolingSchemes{poolingSchemeIndex}.coneInputsNum;
        poolingSchemeName = sprintf('%d-coneInput',conesNumInRFcenter);
        inputConeIndicesForAllCombinations = rfCenterPoolingSchemes{poolingSchemeIndex}.inputConeIndices;
        maxInputConeDistance = rfCenterPoolingSchemes{poolingSchemeIndex}.maxInputConeDistance;
        
        % Generate a family of perfect Gaussian profiles with characteristic radii appropriate
        % for a the n-cone RF center
        characteristicRadiusDegs = 0.5*(maxInputConeDistance+mean(coneAperturesDegs(:)))/3 ;
        characteristicRadiiDegsExamined = characteristicRadiusDegs * logspace(log10(0.5), log10(16), 16);
        GaussiansEnsemble = generateGaussiansEnsemble(thePSFsupportDegsHR, characteristicRadiiDegsExamined, [0 0]);
        
        % Fourier Analysis of Gaussian profiles
        [~, GaussianSFtuningEnsemble, spatialFrequencySupportForGaussianSubregions] = ...
            analyzeGaussianSubregionEnsemble(GaussiansEnsemble, thePSFsupportDegsHR);
            
        % Analyze each input cone combination
        inputConeCombinationsNum = size(inputConeIndicesForAllCombinations,1);
        for inputConeCombinationIndex = 1:inputConeCombinationsNum 
           % fprintf('Analyzing %d of %d spatial input schemes for a %d-cone center subregion.\n',  ...
           %     inputConeCombinationIndex, size(inputConeIndicesForAllCombinations,1), rfCenterPoolingSchemes{poolingSchemeIndex}.coneInputsNum);
      
            % Retrieve the indices of cones feeding into the RF center
            inputConeIndices = inputConeIndicesForAllCombinations(inputConeCombinationIndex,:);
            retinalConeImage = generateRetinalConeImage(inputConeIndices, conePosDegs, coneApertureProfile, thePSFsupportDegsHR);
            
            % Compute the visual cone image = retinal cone image * PSF
            visualConeImage = conv2(retinalConeImage, thePSFHR, 'same');
            
            % Fourier analysis
            [retinalSFtuningMajorAxis, visualSFtuningMajorAxis, ...
             retinalSFtuningMinorAxis, visualSFtuningMinorAxis, ...
             retinalSFSpectrum, visualSFSpectrum, spatialFrequencySupport] = analyzeSpectra(...
                            retinalConeImage, visualConeImage, thePSFsupportDegsHR);
            
            % Checks and memory allocation
            if (inputConeCombinationIndex == 1)
                % Ensure that the spatial frequency support matches that of the Guassian profiles
                assert(any(abs(spatialFrequencySupportForGaussianSubregions - spatialFrequencySupport)<0.0001), ...
                    'spatial frequency axes mismatch');
        
                retinalConeImageSpectrum = zeros(inputConeCombinationsNum, size(retinalSFSpectrum,1), size(retinalSFSpectrum,2));
                visualConeImageSpectrum = zeros(inputConeCombinationsNum, size(visualSFSpectrum,1), size(visualSFSpectrum,2));
                retinalSpatialFrequencyTuningsMajorAxis = zeros(inputConeCombinationsNum, length(retinalSFtuningMajorAxis));
                visualSpatialFrequencyTuningsMajorAxis = zeros(inputConeCombinationsNum, length(visualSFtuningMajorAxis));
                retinalSpatialFrequencyTuningsMinorAxis = zeros(inputConeCombinationsNum, length(retinalSFtuningMinorAxis));
                visualSpatialFrequencyTuningsMinorAxis = zeros(inputConeCombinationsNum, length(visualSFtuningMinorAxis));
            end
            
            % Log data in
            retinalConeImageSpectrum(inputConeCombinationIndex,:,:) = retinalSFSpectrum;
            visualConeImageSpectrum(inputConeCombinationIndex,:,:) = visualSFSpectrum;
            retinalSpatialFrequencyTuningsMajorAxis(inputConeCombinationIndex,:) = retinalSFtuningMajorAxis;
            visualSpatialFrequencyTuningsMajorAxis(inputConeCombinationIndex,:) = visualSFtuningMajorAxis;
            retinalSpatialFrequencyTuningsMinorAxis(inputConeCombinationIndex,:) = retinalSFtuningMinorAxis;
            visualSpatialFrequencyTuningsMinorAxis(inputConeCombinationIndex,:) = visualSFtuningMinorAxis;
        end % inputConeCombinationIndex
        
        % Mean over input cone combinations
        retinalConeImageSpectrum = squeeze(mean(retinalConeImageSpectrum, 1));
        visualConeImageSpectrum = squeeze(mean(visualConeImageSpectrum,1));
        retinalSpatialFrequencyTuningMajorAxis = squeeze(mean(retinalSpatialFrequencyTuningsMajorAxis,1));
        retinalSpatialFrequencyTuningMinorAxis = squeeze(mean(retinalSpatialFrequencyTuningsMinorAxis,1));
        visualSpatialFrequencyTuningMajorAxis = squeeze(mean(visualSpatialFrequencyTuningsMajorAxis,1));
        visualSpatialFrequencyTuningMinorAxis = squeeze(mean(visualSpatialFrequencyTuningsMinorAxis,1));

        switch spatialFrequencyAxisToMatch
            case 'major'
                % Match the major spatial RF axis. This corresponding model will have larger RF centers, and therefore lower peak sensitivity
                spatialFrequencyTuningToMatch = visualSpatialFrequencyTuningMajorAxis;
            case 'minor'
                % Match the minor spatial RF axis. This corresponding model will have smaller RF centers, and therefore higher peak sensitivity
                spatialFrequencyTuningToMatch = visualSpatialFrequencyTuningMinorAxis;
            case 'average'
                 % Match the average of the minor & major spatial RF axis.
                spatialFrequencyTuningToMatch = 0.5*(visualSpatialFrequencyTuningMinorAxis+visualSpatialFrequencyTuningMajorAxis);
        end
        
        
        % Determine the Gaussian whose SF tuning best matches the visualSpatialFrequencyTuningMinorAxis
        [matchingCharacteristicRadiusDegs,  matchingPeakSensitivity, ...
         matchingGaussian, matchingSFrange, ~] = determineMatchingGaussianInFrequencyDomain(characteristicRadiiDegsExamined, ...
            GaussianSFtuningEnsemble, spatialFrequencyTuningToMatch, spatialFrequencySupport, sensitivityRangeOverWhichToMatchSFtuning, ...
            thePSFsupportDegsHR, visualConeImage);
        
        if (visualizeFits)
            % Fourier Analysis of matching Gaussian
            [~, matchingGaussianSFtuning, ~] = analyzeGaussianSubregionEnsemble(...
                reshape(matchingGaussian, [1 size(matchingGaussian,1), size(matchingGaussian,2)]), ...
                thePSFsupportDegsHR);
        
            figNo = poolingSchemeIndex;
            
            % Visualized SF range (c/deg)
            sfRangeVisualized = [0.5 300];
            
            visualizeAnalysis(figNo, quadrantName, subjectID, patchEccRadiusDegs, conesNumInRFcenter, ...
                spatialFrequencySupport, sfRangeVisualized,...
                thePSFsupportDegsHR, thePSFHR, retinalConeImage, visualConeImage, ...
                retinalConeImageSpectrum, visualConeImageSpectrum, ...
                retinalSpatialFrequencyTuningMajorAxis, retinalSpatialFrequencyTuningMinorAxis, ...
                visualSpatialFrequencyTuningMajorAxis, visualSpatialFrequencyTuningMinorAxis, ...
                GaussianSFtuningEnsemble, spatialFrequencySupportForGaussianSubregions, ...
                matchingGaussian, matchingGaussianSFtuning, matchingSFrange,  ...
                overlayMatchingGaussianProfileOnConeImages, exportFig);
        end
        
        % Log data in
        deconvolutionStruct.data(poolingSchemeName) = struct(...
            'characteristicRadiusDegs', matchingCharacteristicRadiusDegs, ...
            'peakSensitivity',  matchingPeakSensitivity ...
        );
    end % poolingSchemeIndex
end


function [matchingCharacteristicRadiusDegs,  matchingPeakSensitivity, ...
         matchingGaussian, matchingSFrange, centerPixelCoords] = determineMatchingGaussianInFrequencyDomain(characteristicRadiiDegsExamined, ...
            GaussianSFtuningEnsemble, targetSFtuning, spatialFrequencySupport, sensitivityRangeOverWhichToMatchSFtuning, ...
            thePSFsupportDegsHR, visualConeImage)

    % Determine the characteristic radius of the Gaussian whose SF
    % tuning best matches the targetSFtuning
    [matchingCharacteristicRadiusDegs, matchingSFrange] = determineMatchingCharacteristicRadius(characteristicRadiiDegsExamined, ...
            GaussianSFtuningEnsemble, targetSFtuning, spatialFrequencySupport, sensitivityRangeOverWhichToMatchSFtuning);
        
    % Generate the matching Gaussian
    [matchingGaussian, matchingPeakSensitivity, centerPixelCoords] = ...
        generateMatchingGaussian(matchingCharacteristicRadiusDegs, thePSFsupportDegsHR, visualConeImage);
end

function [matchingCharacteristicRadiusDegs,  matchingSFrange] = determineMatchingCharacteristicRadius(characteristicRadiiDegsExamined, ...
            GaussianSFtuningEnsemble, targetSFtuning, spatialFrequencySupport, sensitivityRange)
        
    % Find SF range over which to match SF tuning
    idx = find(spatialFrequencySupport>=0);
    GaussianSFtuningEnsemble = GaussianSFtuningEnsemble(:, idx);
    targetSFtuning = targetSFtuning(idx);
    spatialFrequencySupport = spatialFrequencySupport(idx);
    
    % Find index of spatial frequency at which we are closese to the max
    % of the selected sensitivity range. This will be at the low-end of the
    % spatial frequency support
    upperSensitivity = max(sensitivityRange)*max(targetSFtuning);
    [~,idxLow] = min(abs(targetSFtuning-upperSensitivity));
    matchingSFrange(1) = spatialFrequencySupport(idxLow);
    
    % Find index of spatial frequency at which we are closese to the min
    % of the selected sensitivity range. This will be at the high-end of the
    % spatial frequency support. We do not look beyond to the first zero crossing
    lowerSensitivity = min(sensitivityRange)*max(targetSFtuning);
    minDist = Inf;
    idx = idxLow;
    keepGoing = true;
    while (idx <= numel(targetSFtuning)) && (keepGoing)
        dist = abs(targetSFtuning(idx)-lowerSensitivity);
        if (dist < minDist)
            minDist = dist;
            idxHigh = idx;
        end
        if (targetSFtuning(idx) < 0.02*max(targetSFtuning))
            keepGoing = false;
        end
        idx = idx + 1;
    end    
    matchingSFrange(2) = spatialFrequencySupport(idxHigh);

    logSubSamplingSampling = false;
    if (logSubSamplingSampling)
        sampledSpatialFrequencySupport = logspace(log10(matchingSFrange(1)), log10(matchingSFrange(2)), 50);
        sampledSpatialFrequencySupport = sort(unique(sampledSpatialFrequencySupport), 'ascend');

        sfIndicesOfInterest = zeros(1,numel(sampledSpatialFrequencySupport));
        for k = 1:numel(sampledSpatialFrequencySupport)
            targetSF = sampledSpatialFrequencySupport(k);
            [~,sfIndicesOfInterest(k)] = min(abs(spatialFrequencySupport-targetSF));
        end
        sfIndicesOfInterest = sort(unique(sfIndicesOfInterest));
    else
        sfIndicesOfInterest = idxLow:idxHigh;
    end
    
    
    
    GaussianSFtuningEnsemble = GaussianSFtuningEnsemble(:,sfIndicesOfInterest);
    targetSFtuning = targetSFtuning(sfIndicesOfInterest);
    
    m = mean(GaussianSFtuningEnsemble,2);
    [m,idx] = sort(m, 'ascend');
    characteristicRadiiDegsExamined = characteristicRadiiDegsExamined(idx);
    GaussianSFtuningEnsemble = GaussianSFtuningEnsemble(idx,:);
    
    mTarget = mean(targetSFtuning);
    if (mTarget < m(1))
        matchingCharacteristicRadiusDegs = characteristicRadiiDegsExamined(1);
        fprintf(2,'Examined characteristic radii are all too small for this PSF');
        return;
    end
    
    if (mTarget > m(end))
        matchingCharacteristicRadiusDegs = characteristicRadiiDegsExamined(end);
         fprintf(2,'Examined characteristic radii are all too large for this PSF');
        return;
    end
    
    % Find the 2 Gaussian SF tunings that are closest to the target SF tuning
    rmsErrors = sum((bsxfun(@minus, GaussianSFtuningEnsemble, targetSFtuning)).^2,2);
    [rmsErrors,rmsIndices] = sort(rmsErrors);
    characteristicRadiiDegsExamined = characteristicRadiiDegsExamined(rmsIndices);

    % Compute the matchingCharacteristicRadius by weighted sum of the
    % closest 2 characteristic radii
    interpolationIndices = 1:2;
    interpolationWeights(1) = rmsErrors(interpolationIndices(2)) / sum(rmsErrors(interpolationIndices));
    interpolationWeights(2) = rmsErrors(interpolationIndices(1)) / sum(rmsErrors(interpolationIndices));
    matchingCharacteristicRadiusDegs = sum(interpolationWeights .* characteristicRadiiDegsExamined(interpolationIndices));
end

function [matchingGaussian, peakSensitivity, centerPixelCoord] = generateMatchingGaussian(...
            matchingCharacteristicRadiusDegs,thePSFsupportDegsHR, visualConeImage)
        
   % Generate a binary version of the visual cone image (to find the centroid)
   binaryVisualConeImage = visualConeImage/max(visualConeImage(:));
   idx = find(binaryVisualConeImage > 0.01);
   binaryVisualConeImage = 0 * binaryVisualConeImage;
   binaryVisualConeImage(idx) = 1.0;
   
   % Find the centroid
   s = regionprops(binaryVisualConeImage,visualConeImage, 'WeightedCentroid');
   centroids = cat(1,s.WeightedCentroid);
   center = median(centroids,1);
   centerPixelCoord = round(center);
   
   % Generate a Gaussian with the matching characteristic radius at the
   % centroid of the visual cone image
   matchingGaussian = generateGaussiansEnsemble(...
       thePSFsupportDegsHR, matchingCharacteristicRadiusDegs, thePSFsupportDegsHR(centerPixelCoord));
   
   % Make the area of the Gaussian equal to the area of visualConeImage
   gain = sum(visualConeImage(:));
   matchingGaussian = gain * squeeze(matchingGaussian(1,:,:));
   peakSensitivity = max(matchingGaussian(:));
end
        

function [rfSpectra, spatialFrequencyTunings, spatialFrequencySupport] = analyzeGaussianSubregionEnsemble(RFprofiles, thePSFsupportDegsHR)

    % Compute spatial frequency support
    profilesNum = size(RFprofiles,1);
    N = 2^nextpow2(size(RFprofiles,2));
    spatialSampleSize = thePSFsupportDegsHR(2)-thePSFsupportDegsHR(1);
    spatialFrequencySupport = computeSpatialFrequencySupport(spatialSampleSize, N);
    idx = find(spatialFrequencySupport>=0);

    spatialFrequencyTunings = zeros(profilesNum, N);
    rfSpectra = zeros(profilesNum, N, N);
    % Compute 2D and 1D spectra
    for k = 1:profilesNum
        rf = squeeze(RFprofiles(k,:,:));
        rfSpectrum = abs(fftshift(fft2(rf, N,N))/N);
        rfSpectrum = rfSpectrum / max(rfSpectrum(:));
        spatialFrequencyTunings(k,:) = squeeze(rfSpectrum(idx(1),:));
        rfSpectra(k,:,:) = rfSpectrum;
    end
end

function [retinalSpatialFrequencyTuningMajorAxis, visualSpatialFrequencyTuningMajorAxis, ...
          retinalSpatialFrequencyTuningMinorAxis, visualSpatialFrequencyTuningMinorAxis, ...
    retinalConeImageSpectrum, visualConeImageSpectrum, spatialFrequencySupport] = analyzeSpectra(retinalConeImage, visualConeImage, thePSFsupportDegsHR)
    
    % Rotate images so the major axis (longest) is along the x-axis
    [imRotation, ellipseAxesLengthRatio] = determineImageRotation(visualConeImage);
    retinalConeImage = imrotate(retinalConeImage, -(imRotation+90), 'bilinear', 'crop');
    visualConeImage  = imrotate(visualConeImage,  -(imRotation+90), 'bilinear', 'crop');

    % Perform FFT
    N = 2^nextpow2(size(retinalConeImage,2));
    retinalConeImageSpectrum = abs(fftshift(fft2(retinalConeImage, N,N))/N);
    visualConeImageSpectrum  = abs(fftshift(fft2(visualConeImage, N,N))/N);

    % Normalize to max retinal
    maxRetinal = max(retinalConeImageSpectrum(:));
    retinalConeImageSpectrum = retinalConeImageSpectrum / maxRetinal;
    visualConeImageSpectrum = visualConeImageSpectrum  / maxRetinal;
    
    % Compute spatial frequency support
    spatialSampleSize = thePSFsupportDegsHR(2)-thePSFsupportDegsHR(1);
    spatialFrequencySupport = computeSpatialFrequencySupport(spatialSampleSize, N);

    % Retrieve 1D {retinal, visual} slices through the major&minor axes
    idx = find(spatialFrequencySupport>=0);
    retinalSpatialFrequencyTuningMinorAxis = squeeze(retinalConeImageSpectrum(idx(1),:));
    visualSpatialFrequencyTuningMinorAxis = squeeze(visualConeImageSpectrum(idx(1),:));
    retinalSpatialFrequencyTuningMajorAxis = squeeze(retinalConeImageSpectrum(:,idx(1)));
    visualSpatialFrequencyTuningMajorAxis = squeeze(visualConeImageSpectrum(:,idx(1)));
    
%     if (ellipseAxesLengthRatio < 1.1)
%         fprintf('Ellipse axes length ratio < 1.1. Will average minor & major SF tuning functions.\n');
%         % Average the major and minor sf tuning if the ellipse is not elongated
%         retinalSpatialFrequencyTuningMajorAxis = 0.5*(...
%             retinalSpatialFrequencyTuningMinorAxis(:) + retinalSpatialFrequencyTuningMajorAxis(:) ...
%             );
%         visualSpatialFrequencyTuningMajorAxis = 0.5*(...
%             visualSpatialFrequencyTuningMinorAxis(:) + visualSpatialFrequencyTuningMajorAxis(:) ...
%             );
%         retinalSpatialFrequencyTuningMinorAxis = retinalSpatialFrequencyTuningMajorAxis';
%         visualSpatialFrequencyTuningMinorAxis = visualSpatialFrequencyTuningMajorAxis';
%     end
    
end

function spatialFrequencySupport = computeSpatialFrequencySupport(spatialSampleSize, N)
    nyquistSF = 1/(2*spatialSampleSize);
    deltaSF = nyquistSF/(N/2);
    spatialFrequencySupport = linspace(-nyquistSF+deltaSF,nyquistSF, N);
end


function [imRotation, axisLengthRatio] = determineImageRotation(visualConeImage)
    binaryVisualConeImage = visualConeImage*0;
    idx = find(visualConeImage > 0.4*max(visualConeImage(:)));
    binaryVisualConeImage(idx) = 1;
    stats = regionprops(binaryVisualConeImage, 'all');
    axisLengthRatio = stats.MajorAxisLength/stats.MinorAxisLength;
    imRotation = stats.Orientation;        
end


function RFcenterProfiles = generateGaussiansEnsemble(thePSFsupportDegsHR, characteristicRadiiDegs, center)
    
    [X,Y] = meshgrid(thePSFsupportDegsHR,thePSFsupportDegsHR);
    RFcenterProfiles = zeros(numel(characteristicRadiiDegs), size(X,1), size(X,2));
    for k = 1:numel(characteristicRadiiDegs)
        sigma = characteristicRadiiDegs(k);
        GaussianProfile = exp(-((X-center(1))/sigma).^2) .* exp(-((Y-center(2))/sigma).^2);
        % Unit volume
        RFcenterProfiles(k,:,:) = GaussianProfile / sum(GaussianProfile(:));
    end
end

function retinalConeImage = generateRetinalConeImage(inputConeIndices, conePosDegs, coneApertureProfile, thePSFsupportDegs)

    % A 2D grid of delta functions placed at the centers of cones
    retinalConeImage = zeros(length(thePSFsupportDegs), length(thePSFsupportDegs));
    for iCone = 1:numel(inputConeIndices)
         coneXpos = conePosDegs(inputConeIndices(iCone),1);
         coneYpos = conePosDegs(inputConeIndices(iCone),2);
         [~,ix] = min(abs(coneXpos - thePSFsupportDegs));
         [~,iy] = min(abs(coneYpos - thePSFsupportDegs));
         if ( (coneXpos>=thePSFsupportDegs(1)) && (coneXpos<=thePSFsupportDegs(end)) && ...
              (coneYpos>=thePSFsupportDegs(1)) && (coneYpos<=thePSFsupportDegs(end)) )
               retinalConeImage(iy,ix) = 1;
         end
    end
    retinalConeImage = conv2(retinalConeImage, coneApertureProfile, 'same');
end


function poolingSchemes = generateCenterPoolingSchemes(conesNumInRFcenterTested, coneAperturesDegs, conePosDegs)

    poolingSchemes = cell(1, numel(conesNumInRFcenterTested));
    coneCharacteristicRadiusDegs = 0.5*mean(coneAperturesDegs)/3;
    
    for schemeIndex = 1:numel(conesNumInRFcenterTested)
        coneInputsNum = conesNumInRFcenterTested(schemeIndex);
        switch (coneInputsNum)
            case 1
                scenariosNum = 1;
                startingAngle = 0; deltaAngle = 0;
                offsetRadius = 0;
            case 2
                scenariosNum = 6;
                startingAngle = 30; deltaAngle = 60;
                offsetRadius = coneCharacteristicRadiusDegs/WatsonRGCModel.coneApertureToDiameterRatio*3;
            case 3
                scenariosNum = 6;
                startingAngle = 0; deltaAngle = 60;
                offsetRadius = 2/sqrt(3.0)*coneCharacteristicRadiusDegs/WatsonRGCModel.coneApertureToDiameterRatio*3;
            case 4
                scenariosNum = 6;
                startingAngle = 0; deltaAngle = 50;
                offsetRadius = 2/sqrt(3.0)*coneCharacteristicRadiusDegs/WatsonRGCModel.coneApertureToDiameterRatio*3;
            case 5
                scenariosNum = 2;
                startingAngle = 0; deltaAngle = 180;
                offsetRadius = coneCharacteristicRadiusDegs/WatsonRGCModel.coneApertureToDiameterRatio*3;
            otherwise
                scenariosNum = 1;
                startingAngle = 0; deltaAngle = 0;
                offsetRadius = 0;
        end
        
        poolingCenters = zeros(scenariosNum,2);
        inputConeIndices = zeros(scenariosNum,coneInputsNum);
        for iAngle = 1:scenariosNum
            poolingCenters(iAngle,1) = offsetRadius*cosd(startingAngle+iAngle*deltaAngle);
            poolingCenters(iAngle,2) = offsetRadius*sind(startingAngle+iAngle*deltaAngle);
            [~,inputConeIndices(iAngle,:)] = pdist2(conePosDegs, poolingCenters(iAngle,:), 'euclidean', 'smallest', coneInputsNum);
        end
        poolingSchemes{schemeIndex} = struct;
        poolingSchemes{schemeIndex}.coneInputsNum = coneInputsNum; 
        poolingSchemes{schemeIndex}.poolingCenters = poolingCenters;
        poolingSchemes{schemeIndex}.inputConeIndices = inputConeIndices;
        if (coneInputsNum == 1)
            poolingSchemes{schemeIndex}.maxInputConeDistance = mean(coneAperturesDegs);
        else
            poolingSchemes{schemeIndex}.maxInputConeDistance = max(pdist(conePosDegs(inputConeIndices,:)));
        end
    end
end



% ======= ONLY VISUALIZATION ROUTINES BELOW ======

function visualizeAnalysis(figNo, quadrantName, PolansSubjectID, patchEccRadiusDegs, conesNumInRFcenter, spatialFrequencySupport, sfRange, ...
                thePSFsupportDegsHR, thePSFHR, retinalConeImage, visualConeImage, ...
                retinalConeImageSpectrum, visualConeImageSpectrum, ...
                retinalSpatialFrequencyTuningMajorAxis, retinalSpatialFrequencyTuningMinorAxis, ...
                visualSpatialFrequencyTuningMajorAxis, visualSpatialFrequencyTuningMinorAxis, ...
                GaussianSFtuning, spatialFrequencySupportForGaussianSubregions, ...
                matchingGaussian, matchingGaussianSFtuning, matchingSFrange, ...
                overlayMatchingGaussianProfileOnConeImages, exportFig)
            
    xyRange = 0.05 * [-1 1];
    xyTicks = -0.5 : 0.02 : 0.5;
    sfTicks = [0.01 0.03 0.1 0.3 1 3 10 30 100 300];

    redsLUT = brewermap(1024, '*reds');
    bluesLUT = brewermap(1024, '*blues');
    spectralLUT = brewermap(1024, '*greys');
        
    figWidthInches = 20;
    figHeightInches = 18;
    plotlabOBJ = setupPlotLab(0, figWidthInches, figHeightInches);
    
    hFig = figure(figNo); clf;
    set(hFig, 'Name', sprintf('subject %d, %s quadrant, ecc: %2.1f degs', PolansSubjectID, quadrantName, patchEccRadiusDegs));
    rowsNum = 3; colsNum = 3;
    theAxesGrid = plotlab.axesGrid(hFig, ...
            'rowsNum', rowsNum, ...
            'colsNum', colsNum, ...
            'leftMargin', 0.05, ...
            'rightMargin', 0.00, ...
            'widthMargin', 0.08, ...
            'heightMargin', 0.09, ...
            'bottomMargin', 0.05, ...
            'topMargin', 0.01);
    
    % Hide all axes
    for r = 1:rowsNum
        for c = 1:colsNum
            plotlab.toggleAxesVisibility(theAxesGrid{r,c});
        end
    end
    
    % The retinal cone image
    ax = theAxesGrid{1,1}; showYLabel = true; normalizeForVisibility = true;
    renderConeImage(ax, thePSFsupportDegsHR, retinalConeImage, spectralLUT, xyRange, xyTicks, ...
        'retinal cone image', 'r', showYLabel, normalizeForVisibility);
    if (overlayMatchingGaussianProfileOnConeImages)
        hold(ax, 'on');
        overlayContourImage(ax, thePSFsupportDegsHR, matchingGaussian, exp(-1)*[1 1], xyRange, xyTicks);
    end
    
    % The computed visual cone image
    ax = theAxesGrid{2,1}; showYLabel = true; normalizeForVisibility = true;
    renderConeImage(ax, thePSFsupportDegsHR, visualConeImage, spectralLUT, xyRange, xyTicks, ...
        sprintf('visual cone image'),  'b', showYLabel, normalizeForVisibility);
    if (overlayMatchingGaussianProfileOnConeImages)
        hold(ax, 'on');
        overlayContourImage(ax, thePSFsupportDegsHR, matchingGaussian, exp(-1)*[1 1], xyRange, xyTicks);
    end
    
    % The point spread function employed
    ax = theAxesGrid{3,1}; showYLabel = true; normalizeForVisibility = true;
    renderConeImage(ax, thePSFsupportDegsHR, thePSFHR, spectralLUT, xyRange, xyTicks, ...
        sprintf('subject #%d PSF\n(%s quadrant, ecc: %2.1f degs)', PolansSubjectID, quadrantName, patchEccRadiusDegs), ...
        'k', showYLabel, normalizeForVisibility);
   
    
    % 2D specta of retinal cone image
    ax = theAxesGrid{1,2}; showYLabel = true;
    renderSpatialFrequency2Dspectra(ax, spatialFrequencySupport, sfRange, sfTicks, retinalConeImageSpectrum, redsLUT, ...
        sprintf('amplitude spectrum'), 'r', showYLabel);
   
    % SF tuning slices along the major axis (axis of elongation so lower SFs)
    ax = theAxesGrid{1,3}; showXLabel = true; showYLabel = true; reverseXYaxes = true;
    renderSpatialFrequency1Dspectra(ax, spatialFrequencySupport, sfRange, sfTicks, ...
        retinalSpatialFrequencyTuningMajorAxis, visualSpatialFrequencyTuningMajorAxis, 'major axis', reverseXYaxes, showXLabel, showYLabel);
  
    % 2D spectra of visual cone image
    ax = theAxesGrid{2,2}; showYLabel = true;
    renderSpatialFrequency2Dspectra(ax, spatialFrequencySupport, sfRange, sfTicks, visualConeImageSpectrum, bluesLUT, ...
        sprintf('amplitude spectrum'), 'b', showYLabel);
    
    % SF tuning slices along the minor axis - highest SFs
    ax = theAxesGrid{2,3}; showXLabel = true; showYLabel = true; reverseXYaxes = false;
    renderSpatialFrequency1Dspectra(ax, spatialFrequencySupport, sfRange, sfTicks, ...
        retinalSpatialFrequencyTuningMinorAxis, visualSpatialFrequencyTuningMinorAxis, 'minor axis', reverseXYaxes, showXLabel, showYLabel);
    
    % The best matching Gaussian profile
    ax = theAxesGrid{3,2};
    showYLabel = true; normalizeForVisibility = false;
    % Since the visualConeImage is plotted normalized, divide by
    % max(visualConeImage(:)) when plotting the matchingGaussian.
    % This is only for better visualization purposes.
    g = 1/max(visualConeImage(:));
    renderConeImage(ax, thePSFsupportDegsHR,  g*matchingGaussian, spectralLUT, xyRange, xyTicks, ...
       sprintf('matching Gaussian profile (%2.0f-%2.0f cpd)',matchingSFrange(1), matchingSFrange(2)), ...
       'b', showYLabel, normalizeForVisibility);
    hold(ax, 'on');
    overlayContourImage(ax, thePSFsupportDegsHR, matchingGaussian, exp(-1)*[1 1], xyRange, xyTicks);
    
    % 1D spectra of Gaussian profiles with different characteristic radii
    ax = theAxesGrid{3,3};
    renderSpatialFrequency1DspectraForGaussianSubregions(ax, spatialFrequencySupportForGaussianSubregions, ...
        [1 100], sfTicks, GaussianSFtuning, ...
        visualSpatialFrequencyTuningMinorAxis, visualSpatialFrequencyTuningMajorAxis, ...
        matchingGaussianSFtuning, matchingSFrange);
    
    
    % Finish rendering
    drawnow;
    
    % Export to PDF
    if (exportFig)
        isetbioPrefs = getpref('isetbio');
        isetbioRootDir = strrep(isetbioPrefs.validationRootDir, 'validation', '');
        exportDir = fullfile(isetbioRootDir, 'calculators/mosaicConnector/exports');
        pdfFileName = sprintf('Deconv_PolansSID_%d_%s_%2.1fdegs_%d-ConeCenter',  ...
            PolansSubjectID, quadrantName, patchEccRadiusDegs, conesNumInRFcenter);
        plotlabOBJ.exportFig(hFig, 'pdf', pdfFileName, exportDir);
        %setupPlotLab(-1);
    end
    
end


function renderConeImage(ax, thePSFsupportDegsHR, coneImage, imageLUT, xyRange, xyTicks, titleLabel, titleColor, addYLabel, normalizeForVisibility)
    if (normalizeForVisibility)
        coneImage = coneImage/max(coneImage(:));
    end
    
    % Show the axes
    plotlab.toggleAxesVisibility(ax);
    imagesc(ax,thePSFsupportDegsHR, thePSFsupportDegsHR, coneImage); hold on;
    colormap(ax, imageLUT);
    axis(ax, 'xy'); axis(ax,'square'); 
    set(ax, 'XLim', xyRange, 'YLim', xyRange, 'XTick',xyTicks, 'YTick', xyTicks, 'CLim', [0 1]);
    xlabel(ax,'space (degs)');
    if (addYLabel)
        ylabel(ax,'space (degs)');
    end
    title(ax,titleLabel, 'Color', titleColor);
end

function overlayContourImage(ax, thePSFsupportDegsHR, matchingGaussian, contourLevels, xyRange, xyTicks)
    [X,Y] = meshgrid(thePSFsupportDegsHR, thePSFsupportDegsHR);
    normGaussian = matchingGaussian/max(matchingGaussian(:));
    contour(ax, X,Y, normGaussian, contourLevels, 'LineColor', '[0.3 0.3 0]', 'LineWidth', 6);
    contour(ax, X,Y, normGaussian, contourLevels, 'LineColor', [1 1 0], 'LineWidth', 3);
    axis(ax, 'xy'); axis(ax,'square'); 
    set(ax, 'XLim', xyRange, 'YLim', xyRange, 'XTick',xyTicks, 'YTick', xyTicks, 'CLim', [0 1]);
end

function renderSpatialFrequency2Dspectra(ax, spatialFrequencySupport, sfRange, sfTicks, imageSpectrum, imageLUT, ...
    titleLabel, titleColor, addYLabel)

    % Show the axes
    plotlab.toggleAxesVisibility(ax);
    idx = find(spatialFrequencySupport>0);
    spatialFrequencySupportLog = spatialFrequencySupport(idx);
    imageSpectrum = imageSpectrum(idx, idx);
    imagesc(ax, spatialFrequencySupportLog,spatialFrequencySupportLog, imageSpectrum);
    hold(ax, 'on');
    contourf(ax, spatialFrequencySupportLog,spatialFrequencySupportLog, imageSpectrum, ...
        0.05:0.1:1, 'LineColor', [0.3 0.3 0.3]);
    colormap(ax, imageLUT);
    axis(ax, 'square'); axis(ax, 'xy');
    set(ax, 'CLim', [0 1], 'ZLim', [0 1], 'Color', imageLUT(1,:), 'XScale', 'log', 'YScale', 'log', ...
        'XLim', [spatialFrequencySupportLog(1) sfRange(2)], ...
        'YLim', [spatialFrequencySupportLog(1) sfRange(2)], ...
        'XTick', sfTicks, 'YTick', sfTicks);
    xlabel(ax,sprintf('spatial frequency (c/deg)\n (minor axis)'));
    if (addYLabel)
        ylabel(ax,sprintf('spatial frequency (c/deg)\n (major axis)'));
    end
    title(ax,titleLabel, 'Color', titleColor);
end

function renderSpatialFrequency1Dspectra(ax, spatialFrequencySupport, sfRange, sfTicks, ...
    retinalSpatialFrequencyTuning, visualSpatialFrequencyTuning, axisLabel, reverseXYaxes, addXLabel, addYLabel)
    
     % Show the axes
    plotlab.toggleAxesVisibility(ax);
    
    [~,zeroSFindex] = min(abs(spatialFrequencySupport));
    inputConeConfigsNum = size(retinalSpatialFrequencyTuning,1);
    
    if (contains(axisLabel, 'major'))
        lineStyle = '--';
    else
        lineStyle = '-';
    end
    if (reverseXYaxes)
        hold(ax, 'on');
        line(ax, retinalSpatialFrequencyTuning, spatialFrequencySupport, 'Color', 'r', 'LineStyle', '-', 'LineWidth', 4.0); 
        line(ax, visualSpatialFrequencyTuning, spatialFrequencySupport, 'Color', 'b','LineStyle', '-', 'LineWidth', 4.0); 
        line(ax, retinalSpatialFrequencyTuning, spatialFrequencySupport, 'Color', [1 0.5 0.5], 'LineStyle', lineStyle, 'LineWidth', 2.0); 
        line(ax, visualSpatialFrequencyTuning, spatialFrequencySupport, 'Color', 'c','LineStyle', lineStyle, 'LineWidth', 2.0); 
        scatter(ax, (retinalSpatialFrequencyTuning(:,zeroSFindex))', sfRange(1)*ones(1,inputConeConfigsNum), 's',  'MarkerFaceColor', 'r');
        scatter(ax, (visualSpatialFrequencyTuning(:,zeroSFindex))', sfRange(1)*ones(1,inputConeConfigsNum), 's',  'MarkerFaceColor', 'b');
    else
        hold(ax, 'on');
        line(ax, spatialFrequencySupport, retinalSpatialFrequencyTuning, 'Color', 'r', 'LineStyle', '-', 'LineWidth', 4.0); 
        line(ax, spatialFrequencySupport, visualSpatialFrequencyTuning, 'Color', 'b', 'LineStyle', '-', 'LineWidth', 4.0); 
        line(ax, spatialFrequencySupport, retinalSpatialFrequencyTuning, 'Color', [1 0.5 0.5], 'LineStyle', lineStyle, 'LineWidth', 2.0); 
        line(ax, spatialFrequencySupport, visualSpatialFrequencyTuning, 'Color', 'c', 'LineStyle', lineStyle, 'LineWidth', 2.0); 
        scatter(ax,sfRange(1)*ones(1,inputConeConfigsNum), (retinalSpatialFrequencyTuning(:,zeroSFindex))', 's',  'MarkerFaceColor', 'r');
        scatter(ax,sfRange(1)*ones(1,inputConeConfigsNum), (visualSpatialFrequencyTuning(:,zeroSFindex))', 's',  'MarkerFaceColor', 'b');
    end
    
    axis(ax,'square');
    if (reverseXYaxes)
        set(ax, 'YLim', [sfRange(1) sfRange(2)], 'YScale', 'log', 'YTick', sfTicks, ...
            'XLim', [0 1.01], 'XTick', 0:0.2:1.0);
        if (addXLabel)
            ylabel(ax, sprintf('spatial frequency (c/deg)\n (%s)', axisLabel));
        end
        if (addYLabel)
            xlabel(ax, 'amplitude spectrum');
        end
    else 
        set(ax, 'XLim', [sfRange(1) sfRange(2)], 'XScale', 'log', 'XTick', sfTicks, 'YLim', [0 1.01], 'YTick', 0:0.2:1);
        if (addXLabel)
            xlabel(ax, sprintf('spatial frequency (c/deg)\n (%s)', axisLabel));
        end
        if (addYLabel)
            ylabel(ax, 'amplitude spectrum');
        end
    end
end

function renderSpatialFrequency1DspectraForGaussianSubregions(ax, spatialFrequencySupport, sfRange, sfTicks, ...
        GaussianProfileSpatialFrequencyTuning, visualSpatialFrequencyTuningMinorAxis, visualSpatialFrequencyTuningMajorAxis, ...
        matchingGaussianSpatialFrequencyTuning, matchingSFrange)
    
    % Show the axes
    plotlab.toggleAxesVisibility(ax);
    
    % Line colors for each Gaussian profile
    profilesNum = size(GaussianProfileSpatialFrequencyTuning,1);
    lineColors = brewermap(profilesNum*2, 'greys');
    
    % Find the index of SF = 0
    [~,zeroSFindex] = min(abs(spatialFrequencySupport));

    % Plot SF tunings of all profiles
    for k = 1:profilesNum
        line(ax, spatialFrequencySupport, GaussianProfileSpatialFrequencyTuning(k,:), 'Color', lineColors(k,:)/2, 'LineWidth', 3); hold(ax, 'on');
        line(ax, spatialFrequencySupport, GaussianProfileSpatialFrequencyTuning(k,:), 'Color', lineColors(k,:), 'LineWidth', 1);
        scatter(ax,sfRange(1), (GaussianProfileSpatialFrequencyTuning(k,zeroSFindex))', 's',  'MarkerFaceColor', lineColors(k,:));
    end
    
    % Outline the SF tuning of the matching Gaussian in yellow
    line(ax, spatialFrequencySupport, matchingGaussianSpatialFrequencyTuning, 'Color', [0.3 0.3 0], 'LineWidth', 6);
    line(ax, spatialFrequencySupport, matchingGaussianSpatialFrequencyTuning, 'Color', [1 1 0], 'LineWidth', 3);

    % Plot the SF tunings along the minor and major axis
    line(ax, spatialFrequencySupport, visualSpatialFrequencyTuningMinorAxis, 'Color', 'b', 'LineStyle', '-', 'LineWidth', 4);
    line(ax, spatialFrequencySupport, visualSpatialFrequencyTuningMajorAxis, 'Color', 'b', 'LineStyle', '-', 'LineWidth', 4);
    line(ax, spatialFrequencySupport, visualSpatialFrequencyTuningMinorAxis, 'Color', 'c', 'LineStyle', '-', 'LineWidth', 2);
    line(ax, spatialFrequencySupport, visualSpatialFrequencyTuningMajorAxis, 'Color', 'c', 'LineStyle', '--', 'LineWidth', 2);
    
    % Outline the matched SF range
    outline.x = [matchingSFrange(1) matchingSFrange(1) matchingSFrange(2) matchingSFrange(2)];
    outline.y = [0 1 1 0];
    faces = [1 2 3 4];
    vertices = [outline.x(:) outline.y(:)];
    patch(ax, 'Faces', faces, 'Vertices', vertices, 'FaceColor', [1 0.8 0.8], 'EdgeColor', [1 0 0], 'FaceAlpha', 0.2, 'EdgeAlpha', 0.4);
    
     % Plot the SFtuning points at the edge of the analysis range
    [~,sfIndex] = min(abs(spatialFrequencySupport-matchingSFrange(1)));
    plot(spatialFrequencySupport(sfIndex), visualSpatialFrequencyTuningMinorAxis(sfIndex), 'bo', ...
        'MarkerFaceColor', [0.4 1 1], 'MarkerEdgeColor', [0.1 0.8 0.8], 'MarkerSize', 12);
    [~,sfIndex] = min(abs(spatialFrequencySupport-matchingSFrange(2)));
    plot(spatialFrequencySupport(sfIndex), visualSpatialFrequencyTuningMinorAxis(sfIndex), 'bo', ...
        'MarkerFaceColor', [0.4 1 1], 'MarkerEdgeColor', [0.1 0.8 0.8], 'MarkerSize', 12);
    
    
    % Finish plot
    axis(ax,'square');
    set(ax, 'XLim', [sfRange(1) sfRange(2)], 'XScale', 'log', 'XTick', sfTicks, 'YLim', [0 1.01]);
    xlabel(ax, sprintf('spatial frequency (c/deg)'));
    ylabel(ax, 'amplitude spectrum');
end

function plotlabOBJ = setupPlotLab(mode, figWidthInches, figHeightInches)
    if (mode == 0)
        plotlabOBJ = plotlab();
        plotlabOBJ.applyRecipe(...
                'colorOrder', [1 0 0; 0 0 1], ...
                'axesBox', 'off', ...
                'axesTickDir', 'in', ...
                'renderer', 'painters', ...
                'lineMarkerSize', 8, ...
                'axesTickLength', [0.01 0.01], ...
                'legendLocation', 'SouthWest', ...
                'axesFontSize', 18, ...
                'axesFontAngle', 'italic', ...
                'axesLabelFontSizeMultiplier', 1.1, ...
                'figureWidthInches', figWidthInches, ...
                'figureHeightInches', figHeightInches);
    else
        pause(2.0);
        plotlab.resetAllDefaults();
    end
end 