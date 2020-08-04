function generatePolansOpticsDeconvolutionFiles(obj, deconvolutionOpticsParams, varargin)

    defaultEccTested = [0 0.25 0.5 1 1.5 2.0 2.5 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25];
    
    % Parse input
    p = inputParser;
    p.addParameter('eccTested', defaultEccTested);
    p.addParameter('retinalCharacteristicRadii', [0.1]);
    p.parse(varargin{:});

    % Validate the deconvolutionOpticsParams
    obj.validateDeconvolutionOpticsParams(deconvolutionOpticsParams);
    
    eccTested = -(p.Results.eccTested);
    for eccIndex = 1:numel(eccTested)
        ecc = eccTested(eccIndex);
        doIt(obj.psfDeconvolutionDir, ecc, p.Results.retinalCharacteristicRadii, deconvolutionOpticsParams);
    end
end

function doIt(rootDir, cellEcc, retinalCharacteristicRadii, deconvolutionOpticsParams)

    % Whether to visualize the analysis
    visualizeAnalysis = false;
    subjectsNum = numel(deconvolutionOpticsParams.PolansWavefrontAberrationSubjectIDsToCompute);
    quadrantsNum = numel(deconvolutionOpticsParams.quadrantsToCompute);
    
    subjectIDs = deconvolutionOpticsParams.PolansWavefrontAberrationSubjectIDsToCompute;
    quadrants = deconvolutionOpticsParams.quadrantsToCompute;
    
    % Preallocate memory
    retinalRadius = zeros(numel(retinalCharacteristicRadii), quadrantsNum,  subjectsNum);
    visualRadius = retinalRadius;
    visualGain = visualRadius;
   
    useOLDcode = false;
    
    imposedRefractionErrorDiopters = 0;
    
    for qIndex = 1:quadrantsNum

        switch (quadrants{qIndex})
            case 'horizontal'
                % horizontal meridian (nasal)
                eccXrange = cellEcc(1)*[1 1];
                eccYrange = 0*[1 1];
            case 'superior'
                % vertical meridian (superior)
                eccYrange = cellEcc(1)*[1 1];
                eccXrange = 0*[1 1];
            case 'inferior'
                % vertical meridian (inferior)
                eccYrange = -cellEcc(1)*[1 1];
                eccXrange = 0*[1 1];
            otherwise
                error('Unknown Polans quadrant: ''%s''.', eccQuadrant);
        end
        
        for sIndex = 1:subjectsNum
            % Compute the subject PSF at the desired eccentricity
            deconvolutionOpticsParamsForCondition = struct();...
            deconvolutionOpticsParamsForCondition.PolansWavefrontAberrationSubjectIDsToAverage(1) = ...
                deconvolutionOpticsParams.PolansWavefrontAberrationSubjectIDsToCompute(sIndex);
            deconvolutionOpticsParamsForCondition.quadrantsToAverage{1} = ...
                deconvolutionOpticsParams.quadrantsToCompute{qIndex};
            
            if (useOLDcode)
%                 if (cellEcc <= 15)
%                     wavefrontSpatialSamples = 701;
%                 else
%                     wavefrontSpatialSamples = 1001;
%                 end
% 
%                 pupilDiameterMM = 3.0;
%                 wavelengthsListToCompute = [550];
%                 micronsPerDegree = []; % empty so as to compute for each eccentricity
%                 imposedRefractionErrorDiopters = 0;
% 
%                 [hEcc, vEcc, thePSFs, thePSFsupportDegs] = CronerKaplanRGCModel.psfAtEccentricity(subjectID, ...
%                     imposedRefractionErrorDiopters, pupilDiameterMM, wavelengthsListToCompute, micronsPerDegree, ...
%                     wavefrontSpatialSamples, eccXrange, eccYrange, deltaEcc);
%                 fprintf('PSF computed');
% 
%                 % Make a large PSF (via zero padding to be used to convolve with
%                 % the largest stimuli)
%                 thePSFOriginal = squeeze(thePSFs(1, 1, 1,:,:));
%                 thePSFsupportDegsOriginal = thePSFsupportDegs;
% 
%                 psfSize = size(thePSFOriginal,1);
%                 largePSFsize = 1201;
%                 theLargePSF = zeros(largePSFsize, largePSFsize);
%                 margin = 0.5*(largePSFsize-psfSize);
%                 theLargePSF(margin+(1:psfSize), margin+(1:psfSize)) = thePSFOriginal;
% 
%                 % Make corresponding large PSF support
%                 dS = thePSFsupportDegsOriginal(2)-thePSFsupportDegsOriginal(1);
%                 theLargePSFsupportDegs = -(0.5*(largePSFsize-1)*dS):dS:(0.5*(largePSFsize-1)*dS);
            end
            
            
            % Convolve different retinal pooling regions and compute the
            % visually-mapped pooling region
            
            % parfor retinalRadiusIndex = 1:numel(retinalCharacteristicRadii)
            for retinalRadiusIndex = 1:numel(retinalCharacteristicRadii)
                fprintf('Computing spread for quadrant %s, radius %d/%d\n', ...
                    quadrants{qIndex}, retinalRadiusIndex, numel(retinalCharacteristicRadii));
                
                retinalCharacteristicRadiusDegs = retinalCharacteristicRadii(retinalRadiusIndex);

                [rfSigmasRetinal, rfGainRetinal, rfSigmasVisual, rfGainVisual] = generateConeBasedGaussiaSubregionProfileAndOptics(...
                    eccXrange(1), eccYrange(1), imposedRefractionErrorDiopters, deconvolutionOpticsParamsForCondition, ...
                    retinalCharacteristicRadiusDegs);
                    
                % Saved data
                retinalRadius(retinalRadiusIndex,qIndex, sIndex) = min(rfSigmasRetinal);
                visualRadius(retinalRadiusIndex,qIndex, sIndex) = min(rfSigmasVisual);
                visualGain(retinalRadiusIndex,qIndex, sIndex) = rfGainVisual/rfGainRetinal;
                
                if (visualizeAnalysis)
                    % Extract profile at peak
                    row = round(size(rfPoolingInRetinalSpace,1)/2);
                    rfPoolingRetinalSpaceProfile = squeeze(rfPoolingInRetinalSpace(row,:));
                    rfPoolingVisualSpaceProfile = squeeze(rfPoolingInVisualSpace(row,:));
                    maxX = 0.1;
                    visualizeFit(thePSFsupportDegs, thePSF, subjectID, eccXrange, eccYrange, ...
                        rfPoolingInRetinalSpaceNorm, rfPoolingInVisualSpaceNorm, ...
                        rfPoolingRetinalSpaceProfile, rfPoolingVisualSpaceProfile, ...
                        retinalRadius(retinalRadiusIndex,qIndex, sIndex), visualRadius(retinalRadiusIndex,qIndex, sIndex), ...
                        visualGain(retinalRadiusIndex,qIndex, sIndex), ellipseInRetinalSpace, ellipseInVisualSpace, cMap, maxX);
                end
            end
        end % sIndex
    end % qIndex
    
    % Save data
    dataFileName = fullfile(rootDir, sprintf('ecc_%2.1f_deconvolutions_refractionError_%2.2fD.mat', cellEcc(1), imposedRefractionErrorDiopters));
    fprintf('Saving data to %s\n', dataFileName);
    save(dataFileName, 'retinalPoolingRadii', 'retinalRadius', 'visualRadius', 'visualGain', 'subjectIDs', 'quadrants');   
end

function [rfSigmasRetinal, rfGainRetinal, rfSigmasVisual, rfGainVisual] = ...
    generateConeBasedGaussiaSubregionProfileAndOptics(eccXrange, eccYrange, ...
    imposedRefractionErrorDiopters, deconvolutionOpticsParams, ...
    retinalCharacteristicRadiusDegs)

    % Generate cone positions appropriate for the eccentricity
    patchEccDegs = [eccXrange(1), eccYrange(1)];
    [conePosDegs, coneAperturesDegs, thePSF, thePSFsupportDegs] = ...
        generateConePositionsAndPSFForPatchAtEccentricity(patchEccDegs, ...
        imposedRefractionErrorDiopters, deconvolutionOpticsParams);
    
    
    % Rectangular mesh
    [X,Y] = meshgrid(thePSFsupportDegs, thePSFsupportDegs);
    
    % Cone aperture profile (flat-top)
    coneCharacteristicRadiusDegs = 0.5*mean(coneAperturesDegs)/3;
    coneApertureProfile = exp(-(X/coneCharacteristicRadiusDegs).^2) .* exp(-(Y/coneCharacteristicRadiusDegs).^2);
    flatTopThreshold = exp(-3)*exp(-3);
    coneApertureProfile(coneApertureProfile>=flatTopThreshold) = flatTopThreshold;
    coneApertureProfile(coneApertureProfile<flatTopThreshold) = 0;
    coneApertureProfile = coneApertureProfile / flatTopThreshold;
    
    % 2D cone aperture grid
    deltaFunctionMap2D = zeros(length(thePSFsupportDegs), length(thePSFsupportDegs));
    for iCone = 1:size(conePosDegs,1)
        [~,ix] = min(abs(conePosDegs(iCone,1) - thePSFsupportDegs));
        [~,iy] = min(abs(conePosDegs(iCone,2) - thePSFsupportDegs));
        deltaFunctionMap2D(iy,ix) = 1;
    end
    
    % Cone aperture sampling map
    coneSamplingMap = conv2(deltaFunctionMap2D, coneApertureProfile, 'same');
    
    
    % The continuous activation (perfect Gaussian subregion)
    continuousActivationMapInRetinalSpace = exp(-(X/retinalCharacteristicRadiusDegs).^2).*exp(-(Y/retinalCharacteristicRadiusDegs).^2);
    
    % The continous activation in visual space
    continuousActivationMapInVisualSpace = conv2(continuousActivationMapInRetinalSpace, thePSF, 'same');
    
    % Continuous activation of cone apertures mapped in retinal space
    discreteConeActivationMapInRetinalSpace = coneSamplingMap .* continuousActivationMapInRetinalSpace;
    
    % Continous activation of cone apertures mapped in visual space
    discreteConeActivationMapInVisualSpace =  coneSamplingMap .* continuousActivationMapInVisualSpace;
    
    % Discrete cone activation maps
    integratedConeActivationMapInRetinalSpace = zeros(1, size(conePosDegs,1));
    integratedConeActivationMapInVisualSpace = zeros(1, size(conePosDegs,1));
    
    % Compute cone activations in retinal and visual space by integrating
    % the coneBasedSubregion maps over the cone aperture
    for iCone = 1:size(conePosDegs,1)
        xr = conePosDegs(iCone,1) + coneAperturesDegs(iCone)*0.5*cosd(0:20:360);
        yr = conePosDegs(iCone,2) + coneAperturesDegs(iCone)*0.5*sind(0:20:360);
        [in,on] = inpolygon(X(:), Y(:), xr, yr);
        allIndices = find((in==true) | (on == true));
        if (numel(allIndices) == 0)
            integratedConeActivationMapInRetinalSpace(iCone) = 0;
            integratedConeActivationMapInVisualSpace(iCone) = 0;
        else
            integratedConeActivationMapInRetinalSpace(iCone) = mean(discreteConeActivationMapInRetinalSpace(allIndices));
            integratedConeActivationMapInVisualSpace(iCone) = mean(discreteConeActivationMapInVisualSpace(allIndices));
        end
    end
    
    
    maxActivations = max([max(integratedConeActivationMapInVisualSpace) max(integratedConeActivationMapInRetinalSpace)]);
    integratedConeActivationMapInRetinalSpace = integratedConeActivationMapInRetinalSpace / maxActivations;
    integratedConeActivationMapInVisualSpace = integratedConeActivationMapInVisualSpace/ maxActivations;
    
    

    visualizeSpaceRange = coneCharacteristicRadiusDegs*3*2*10*[-1 1];
    
    figure(15); clf;
    ax = subplot(2,4,1);
    % Render the PSF and the cone apertures
    plotConeApertures(ax, conePosDegs, coneAperturesDegs, zeros(size(coneAperturesDegs))); hold on
    contour(ax, X,Y, thePSF/max(thePSF(:)), 0:0.2:1, 'LineColor', 'r', 'LineWidth', 1.5);
    colormap(ax, gray);
    axis(ax, 'xy');
    axis(ax, 'square');
    set(ax, 'XLim', visualizeSpaceRange, 'YLim', visualizeSpaceRange, 'CLim', [0 1]);
 
    
    
    ax = subplot(2,4,2);
    visualizeContinuousActivationAndPointSpreadFunction(ax, thePSFsupportDegs, visualizeSpaceRange, conePosDegs, coneAperturesDegs,  continuousActivationMapInRetinalSpace, [],  'retinally mapped cone activation')
%     
    % Plot the model subregion (Gaussian) and the PSF
    ax = subplot(2,4,6);
    visualizeContinuousActivationAndPointSpreadFunction(ax, thePSFsupportDegs, visualizeSpaceRange, conePosDegs, coneAperturesDegs, continuousActivationMapInVisualSpace, thePSF, 'visually mapped cone activation')

    
    % Plot the continous cone activation in retinal space
    ax = subplot(2,4,3);
    visualizeDiscreteConeApertureActivationMap(ax, thePSFsupportDegs, visualizeSpaceRange, discreteConeActivationMapInRetinalSpace, 'retinally mapped cone activation')
    
    % Plot the continous cone activation in visual space
    ax = subplot(2,4,7);
    visualizeDiscreteConeApertureActivationMap(ax, thePSFsupportDegs, visualizeSpaceRange, discreteConeActivationMapInVisualSpace, 'visually mapped cone activation');
    drawnow
    
    
    % Fit 2D activation maps
    [fittedRetinalActivationMap, fittedRetinalSlice, hiResPSFsupportDegs, rfSigmasRetinal, rfGainRetinal] = ...
         fitActivationMap(integratedConeActivationMapInRetinalSpace, ...
         conePosDegs(:,1), conePosDegs(:,2), thePSFsupportDegs);
     
    [fittedVisualActivationMap, fittedVisualSlice, ~, rfSigmasVisual, rfGainVisual] = ...
         fitActivationMap(integratedConeActivationMapInVisualSpace, ...
         conePosDegs(:,1), conePosDegs(:,2), thePSFsupportDegs);
     
    
    % Plot the integrated cone activations as they would be measured without optics (retinal view)
    ax = subplot(2,4,4);
    visualizeIntegratedConeApertureActivationMap(ax, visualizeSpaceRange, conePosDegs, coneAperturesDegs, ...
        integratedConeActivationMapInRetinalSpace, fittedRetinalActivationMap, hiResPSFsupportDegs, 'retinally mapped cone activation');
    
    % Plot the integrated cone activations as they would be measured with optics (visual space view)
    ax = subplot(2,4,8);
    visualizeIntegratedConeApertureActivationMap(ax, visualizeSpaceRange, conePosDegs, coneAperturesDegs, ...
        integratedConeActivationMapInVisualSpace,  fittedVisualActivationMap, hiResPSFsupportDegs, 'visually mapped cone activation');
      
    % Render activation of cones (retinal and visual space view) along the horizontal meridian
    ax = subplot(2,4,5);
    midRowCones = find(abs(conePosDegs(:,2)) < 1e-3);
    hold(ax,'on')
    for iCone = 1:numel(midRowCones)
       scatter(ax, conePosDegs(midRowCones,1), integratedConeActivationMapInRetinalSpace(midRowCones), 'o', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', [0 0 0]);
       scatter(ax, conePosDegs(midRowCones,1), integratedConeActivationMapInVisualSpace(midRowCones),  'o', 'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', [1 0 0]);
    end
    plot(ax, hiResPSFsupportDegs , fittedRetinalSlice, 'k-');
    plot(ax, hiResPSFsupportDegs , fittedVisualSlice, 'r-');
    axis(ax, 'square');
    set(ax, 'XLim', visualizeSpaceRange, 'YLim', [0 1]);
    
    drawnow;
    
end

function [fittedActivationMapFull2D, fittedActivationSlice, hiResPSFsupportDegs, rfSigmas, rfGain] = fitActivationMap(...
    activationMap, X, Y, thePSFsupportDegs)

    deltaX = thePSFsupportDegs(2)-thePSFsupportDegs(1);
    hiResPSFsupportDegs = min(thePSFsupportDegs):deltaX/4:max(thePSFsupportDegs);
    
    [rfCenter, rfSigmas, rfGain, rfRotationDegs, rfExponent,rfFunction] = ...
        fitElliptical2DGausianToRF(X(:), Y(:), activationMap(:), deltaX/10, [0 0]);
    
    [XXX,YYY] = meshgrid(hiResPSFsupportDegs, hiResPSFsupportDegs);
    xyData(:,1) = XXX(:);
    xyData(:,2) = YYY(:);
    fittedParams(1) = rfGain;
    fittedParams(2) = rfCenter(1);
    fittedParams(3) = rfCenter(2);
    fittedParams(4) = rfSigmas(1);
    fittedParams(5) = rfSigmas(2);
    fittedParams(6) = rfRotationDegs;
    fittedParams(7) = rfExponent;
    
    fittedActivationMap = rfFunction(fittedParams,xyData);
    fittedActivationMapFull2D = reshape(fittedActivationMap, [numel(hiResPSFsupportDegs) numel(hiResPSFsupportDegs)]);
    
    midRow = find(abs(hiResPSFsupportDegs) < 1e-3);
    fittedActivationSlice = fittedActivationMapFull2D(midRow,:);
end

function  visualizeContinuousActivationAndPointSpreadFunction(ax, thePSFsupportDegs, visualizeSpaceRange, conePosDegs, coneAperturesDegs, continuousActivation, thePSF, titleLabel)
    hold(ax, 'on');
    imagesc(ax, thePSFsupportDegs, thePSFsupportDegs, continuousActivation/max(continuousActivation(:)));
    plotConeApertures(ax, conePosDegs, coneAperturesDegs, zeros(size(coneAperturesDegs)));
    if (~isempty(thePSF))
        [X,Y] = meshgrid(thePSFsupportDegs, thePSFsupportDegs);
        contour(ax, X,Y, thePSF/max(thePSF(:)), 0:0.2:1, 'LineColor', 'r', 'LineWidth', 1.5);
    end
    set(ax, 'XLim', visualizeSpaceRange, 'YLim', visualizeSpaceRange, 'CLim', [0 1]);
    axis(ax, 'xy');
    axis(ax,'square');
    title(ax, sprintf('%s\n(continuous)', titleLabel));
end

function visualizeDiscreteConeApertureActivationMap(ax, thePSFsupportDegs, visualizeSpaceRange, coneBasedSubgegionActivation, titleLabel)
    coneBasedSubgegionActivation = coneBasedSubgegionActivation / max(coneBasedSubgegionActivation(:));
    imagesc(ax, thePSFsupportDegs, thePSFsupportDegs, coneBasedSubgegionActivation); hold on;
    set(ax, 'XLim', visualizeSpaceRange, 'YLim', visualizeSpaceRange);
    axis(ax, 'xy');
    axis(ax, 'square');
    set(ax, 'XLim', visualizeSpaceRange, 'YLim', visualizeSpaceRange, 'CLim', [0 1]);
    title(ax, sprintf('%s\n(discrete)', titleLabel));
end

function visualizeIntegratedConeApertureActivationMap(ax, visualizeSpaceRange, conePosDegs, coneAperturesDegs, ...
    integratedConeBasedSubgegionActivation, fittedActivationMap, hiResPSFsupportDegs, titleLabel)

    hold(ax, 'on');
    %integratedConeBasedSubgegionActivation = integratedConeBasedSubgegionActivation / max(integratedConeBasedSubgegionActivation(:));
    plotConeApertures(ax, conePosDegs, coneAperturesDegs, integratedConeBasedSubgegionActivation);
    [X,Y] = meshgrid(hiResPSFsupportDegs, hiResPSFsupportDegs);
    contour(X,Y, fittedActivationMap, 0.1:0.2:1.0, 'LineColor', 'r', 'LineWidth', 1.5);
    axis(ax, 'xy');
    axis(ax, 'square');
    set(ax, 'Color', [1 1 1], 'XLim', visualizeSpaceRange, 'YLim', visualizeSpaceRange, 'CLim', [0 1], 'ZLim', [0 1]);
    title(ax, sprintf('%s\n(integrated)', titleLabel));
end


function plotConeApertures(ax, conePosDegs, coneAperturesDegs, coneActivations)
    for iCone = 1:size(conePosDegs,1)
        xr = conePosDegs(iCone,1) + coneAperturesDegs(iCone)*0.5*cosd(0:20:360);
        yr = conePosDegs(iCone,2) + coneAperturesDegs(iCone)*0.5*sind(0:20:360);
        patch(ax, xr, yr, [1 1 1]-coneActivations(iCone)*[0 1 1], 'FaceAlpha', 1);
    end
end


function [conePosDegs, coneAperturesDegs, thePSF, thePSFsupportDegs] = ...
    generateConePositionsAndPSFForPatchAtEccentricity(patchEccDegs, imposedRefractionErrorDiopters, deconvolutionOpticsParams)

    recomputeConeMosaic = true;
    recomputeRGCmosaic = false;
    recomputeOptics = true;
    
    patchEccMicrons = 1000*WatsonRGCModel.rhoDegsToMMs(patchEccDegs);

    runParams.rgcMosaicPatchEccMicrons = patchEccMicrons;
    runParams.rgcMosaicPatchSizeMicrons = [0 0];  % Just the max surround region as computed in mosaicsAndOpticsForEccentricity
    runParams.orphanRGCpolicy = 'steal input';
    runParams.maximizeConeSpecificity = 100;
    runParams.noLCA = false;
    runParams.noOptics = false;
    runParams.imposedRefractionErrorDiopters = imposedRefractionErrorDiopters;
    runParams.deconvolutionOpticsParams = deconvolutionOpticsParams;

    [theConeMosaic, ~, theOI] = ...
        mosaicsAndOpticsForEccentricity(runParams, recomputeConeMosaic, recomputeRGCmosaic, recomputeOptics, ...
        '', false);
   
    % Get the cone positions and apertures
    cmStructSerialized = theConeMosaic.geometryStructAlignedWithSerializedConeMosaicResponse();
    conePosDegs = WatsonRGCModel.rhoMMsToDegs(cmStructSerialized.coneLocsMicrons/1000);
    coneAperturesDegs = cmStructSerialized.coneApertures;
    
    % Get PSF slice at target wavelength
    theOptics = oiGet(theOI, 'optics');
    targetWavelength = 550;
    thePSF = opticsGet(theOptics,'psf data',targetWavelength);
    
    % Extract support in arcmin
    psfSupportMicrons = opticsGet(theOptics,'psf support','um');
    if (isfield(theOptics, 'micronsPerDegree'))
        micronsPerDegree = theOptics.micronsPerDegree;
    else
        focalLengthMeters = opticsGet(theOptics, 'focalLength');
        focalLengthMicrons = focalLengthMeters * 1e6;
        micronsPerDegree = focalLengthMicrons * 2 * tand(0.5);
    end

    xGridDegs = psfSupportMicrons{1}/micronsPerDegree;
    thePSFsupportDegs = xGridDegs(1,:);
end


function visualizeFit(thePSFsupportDegs, thePSF, subjectID, eccXrange, eccYrange, ...
                rfPoolingInRetinalSpace, rfPoolingInVisualSpace, ...
                rfPoolingRetinalSpaceProfile, rfPoolingVisualSpaceProfile, ...
                retinalRadius, visualRadius, ...
                visualGain, ellipseInRetinalSpace, ellipseInVisualSpace, cMap, maxX)
            
        % Visualize them
        hFig = figure(1); clf;

        subplotPosVectors = NicePlot.getSubPlotPosVectors(...
               'rowsNum', 2, ...
               'colsNum', 4, ...
               'heightMargin',  0.03, ...
               'widthMargin',    0.02, ...
               'leftMargin',     0.03, ...
               'rightMargin',    0.00, ...
               'bottomMargin',   0.04, ...
               'topMargin',      0.02);

        xLims = maxX*[-1 1];
        % Plot the PSF
        
        ax = subplot('Position', subplotPosVectors(1,1).v);
        theTitle = sprintf('subject #%d PSF at %2.1f, %2.1f degs', subjectID, eccXrange(1), eccYrange(1));
        renderPSF(ax, thePSFsupportDegs, thePSF, xLims, cMap, theTitle)

        ax = subplot('Position', subplotPosVectors(1,2).v);
        renderKernel(ax, thePSFsupportDegs, rfPoolingInRetinalSpace, xLims, cMap, 'pooling in retinal space');
        
        ax = subplot('Position', subplotPosVectors(2,2).v);
        renderKernelProfile(ax, thePSFsupportDegs, rfPoolingRetinalSpaceProfile, [], retinalRadius, [], xLims);
        
        ax = subplot('Position', subplotPosVectors(1,3).v);
        renderKernel(ax, thePSFsupportDegs, rfPoolingInVisualSpace, xLims, cMap, 'pooling in visual space');
            
        ax = subplot('Position', subplotPosVectors(2,3).v);
        renderKernelProfile(ax, thePSFsupportDegs, [], rfPoolingVisualSpaceProfile, [], visualRadius, xLims);

        ax = subplot('Position', subplotPosVectors(1,4).v);
        renderFittedEllipses(ax, ellipseInRetinalSpace, ellipseInVisualSpace, xLims, 'fitted ellipses');
        
        ax = subplot('Position', subplotPosVectors(2,4).v);
        renderKernelProfile(ax, thePSFsupportDegs, rfPoolingRetinalSpaceProfile, rfPoolingVisualSpaceProfile, [], [], xLims);

        
        ax = subplot('Position', subplotPosVectors(2,1).v);
        bar(ax, 0, visualGain, 1, 'r');
        set(ax, 'XLim', [-0.1 0.1], 'YLim', [0 1], 'XTick', [0], 'YTick', 0:0.2:1);
        axis(ax, 'square');
        drawnow;
        
end

function renderPSF(ax, thePSFSupportDegs, thePSF, xLims, cMap, theTitle)
    peakV = max(thePSF(:));
    imagesc(ax, thePSFSupportDegs, thePSFSupportDegs, thePSF); hold(ax, 'on');
    zLevels = logspace(log10(0.05), log10(1), 6)*peakV;
    [X,Y] = meshgrid(thePSFSupportDegs, thePSFSupportDegs);
    contour(ax, X,Y, thePSF, zLevels, 'LineColor', [0.5 0.5 0.5]);
    midRow = round(size(thePSF,1)/2)+1;
    profile = thePSF(midRow,:);
    deltaX = xLims(2)-xLims(1);
    area(ax,thePSFSupportDegs, xLims(1) + deltaX*profile/peakV, xLims(1), 'FaceColor', [1 0.5 0.5], 'EdgeColor', [1 0 0]);
    %plot(ax, [0 0], [-1 1], 'k-'); plot(ax, [-1 1], [0 0], 'k-');
    axis(ax, 'square');  axis(ax, 'xy');
    set(ax, 'XLim', xLims, 'YLim', xLims, 'CLim', [0 peakV]);
    set(ax, 'XTickLabel', {}, 'YTickLabel', {});
    set(ax, 'XTick', -1:0.05:1, 'YTick', -1:0.05:1, 'XTickLabel', {}, 'YTickLabel', {});
    colormap(ax,cMap);
    title(ax, theTitle);
end

function renderKernel(ax, thePSFsupportDegs, rfPoolingInRetinalSpace, xLims, cMap, theTitle)
        contourf(ax, thePSFsupportDegs, thePSFsupportDegs, rfPoolingInRetinalSpace, 6);  hold(ax, 'on');
        plot(ax, [0 0], [-1 1], 'k-'); plot(ax, [-1 1], [0 0], 'k-');
        axis(ax, 'square');  axis(ax, 'xy');
        set(ax, 'XLim', xLims, 'YLim', xLims, 'CLim', [0 1.1]);
        colormap(ax,cMap);
        ticks = -1:0.05:1;
        set(ax, 'XTickLabel', {}, 'YTickLabel', {}, ...
            'XTick', ticks, 'YTick', ticks, ...
            'XTickLabel', ticks*60, 'YTick', ticks*60);
        axis(ax, 'square');
        title(ax, theTitla@Qe);
end


function renderKernelProfile(ax, thePSFsupportDegs, rfPoolingRetinalSpaceProfile, rfPoolingVisualSpaceProfile, retinalRadius, visualRadius, xLims)
        hold(ax, 'on');
        if (~isempty(rfPoolingVisualSpaceProfile))
            area(ax, thePSFsupportDegs, rfPoolingVisualSpaceProfile, 'LineWidth', 1.0, 'FaceColor', [1 0.5 0.5]);
        end
        if (~isempty(rfPoolingRetinalSpaceProfile))
            area(ax, thePSFsupportDegs, rfPoolingRetinalSpaceProfile,  'LineWidth', 1.0, 'FaceColor', [0.5 0.5 0.5]);
        end
        if (~isempty(retinalRadius))
            plot(retinalRadius*[-1 1], max(rfPoolingRetinalSpaceProfile)*exp(-1)*[1 1], 'k-');
        end
        if (~isempty(visualRadius))
            [~,idx] = max(rfPoolingVisualSpaceProfile);
            xo = thePSFsupportDegs(idx);
            plot(visualRadius*[-1 1]+xo, max(rfPoolingVisualSpaceProfile)*exp(-1)*[1 1], 'k-');
        end
        
        plot(ax, [0 0], [0 1], 'k-');
        axis(ax, 'square'); axis(ax, 'xy');
        set(ax, 'XLim', xLims, 'YLim', [0 1]);
        ticks = -1:0.05:1;
        set(ax, 'XTick', ticks,  'XTickLabel', ticks*60, 'YTick', 0:0.2:1, 'YTickLabel', {});
        axis(ax, 'square');
end

function renderFittedEllipses(ax, ellipseInRetinalSpace, ellipseInVisualSpace, xLims, theTitle)
        plot(ax, ellipseInRetinalSpace.x, ellipseInRetinalSpace.y, 'k-', 'LineWidth', 1.5);
        hold(ax, 'on');
        plot(ax, ellipseInVisualSpace.x, ellipseInVisualSpace.y, 'r-', 'LineWidth', 1.5);
        plot(ax, [0 0], [-1 1], 'k-'); plot(ax, [-1 1], [0 0], 'k-');
        axis(ax, 'square');  axis(ax, 'xy');
        ticks = -1:0.05:1;
        set(ax, 'XLim', xLims, 'YLim', xLims,  'XTick', ticks, 'YTick', ticks, ...
            'XTickLabel', ticks*60, 'YTickLabel', ticks*60);
        set(ax, 'XTickLabel', {}, 'YTickLabel', {});
        title(ax, theTitle);
end


function plotlabOBJ = setupPlotLab()
    plotlabOBJ = plotlab();
    plotlabOBJ.applyRecipe(...
            'colorOrder', [0 0 0], ...
            'lineColor', [0.5 0.5 0.5], ...
            'lineWidth', 1.0, ...
            'scatterMarkerEdgeColor', [0.3 0.3 0.9], ...
            'lineMarkerSize', 10, ...
            'axesBox', 'off', ...
            'axesTickDir', 'in', ...
            'renderer', 'painters', ...
            'axesTickLength', [0.01 0.01], ...
            'legendLocation', 'SouthWest', ...
            'figureWidthInches', 24, ...
            'figureHeightInches', 12);
end