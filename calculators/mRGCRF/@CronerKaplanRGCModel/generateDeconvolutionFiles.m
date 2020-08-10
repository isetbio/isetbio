function generateDeconvolutionFiles(obj, deconvolutionOpticsParams, varargin)

    defaultEccTested = [0 0.25 0.5 1 1.5 2.0 2.5 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25];
    
    % Parse input
    p = inputParser;
    p.addParameter('eccTested', defaultEccTested);
    p.parse(varargin{:});

    % Validate the deconvolutionOpticsParams
    obj.validateDeconvolutionOpticsParams(deconvolutionOpticsParams);
    
    eccTested = -(p.Results.eccTested);
    imposedRefractionErrorDiopters = 0;
    
    for eccIndex = 1:numel(eccTested)
        ecc = eccTested(eccIndex);
        doIt(obj, ecc, deconvolutionOpticsParams, imposedRefractionErrorDiopters);
    end
    
end

function doIt(obj, patchEccDegs,  deconvolutionOpticsParams, imposedRefractionErrorDiopters)
    
    % Extra deconvolution optics params
    subjectIDs = deconvolutionOpticsParams.PolansWavefrontAberrationSubjectIDsToCompute;
    quadrants = deconvolutionOpticsParams.quadrantsToCompute;
    subjectsNum = numel(subjectIDs);
    quadrantsNum = numel(quadrants);
   
    useOLDcode = false;
    
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
            patchEccDegs, ...
            'coneMosaicResamplingFactor', 1);
    
        for sIndex = 1:subjectsNum
            % Compute the Polans subject PSF at the desired eccentricity
            PolansSubjectID = deconvolutionOpticsParams.PolansWavefrontAberrationSubjectIDsToCompute(sIndex);
            
            [thePSF, thePSFsupportDegs] = generatePSFForPatchAtEccentricity(PolansSubjectID, ...
                imposedRefractionErrorDiopters, wavelengthSampling, ...
                micronsPerDegree, patchEccDegs);
            
            if (qIndex == 1) && (sIndex == 1)
                obj.setupPlotLab(0, 20, 12);
                resetPlotLabOnExit = true;
            end
            
            if (useOLDcode)
%                 if (patchEcc <= 15)
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
            
            
            % Convolve different retinal pooling regions and compute the visually-mapped pooling region
            
            deconvolutionStruct{ qIndex, sIndex} = analyzeRetinalAndVisualSpatialPoolingSpread(conePosDegs, coneAperturesDegs, ...
                thePSF, thePSFsupportDegs, PolansSubjectID, patchEccDegs);
                
            deconvolutionStruct{qIndex, sIndex}
                
                % Save data

                
                
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
    dataFileName = fullfile(obj.psfDeconvolutionDir, sprintf('ecc_%2.1f_%2.1f_deconvolutions_refractionError_%2.2fD.mat', eccXrange(1), eccYrange(1), imposedRefractionErrorDiopters));
    fprintf('Saving data to %s\n', dataFileName);
    save(dataFileName, ...
            'deconvolutionStruct', ...
            'subjectIDs', 'quadrants');  
        
    if (resetPlotLabOnExit)
        obj.setupPlotLab(-1);
    end
    
end

function deconvolutionStruct = analyzeRetinalAndVisualSpatialPoolingSpread(conePosDegs, ...
    coneAperturesDegs, thePSF, thePSFsupportDegs, ...
    PolansSubjectID, patchEccDegs)
    
    % Rectangular mesh
    [X,Y] = meshgrid(thePSFsupportDegs, thePSFsupportDegs);
    
    % Cone aperture profile (flat-top)
    coneCharacteristicRadiusDegs = 0.5*mean(coneAperturesDegs)/3;
    
    apertureShape = 'disk';
    if (strcmp(apertureShape, 'disk'))
        r = sqrt(X.^2+Y.^2);
        coneApertureProfile = zeros(size(X));
        coneApertureProfile(r<=0.5*mean(coneAperturesDegs)) = 1;
    elseif (strcmp(apertureShape, 'flattopGaussian'))
        coneApertureProfile = exp(-(X/coneCharacteristicRadiusDegs).^2) .* exp(-(Y/coneCharacteristicRadiusDegs).^2);
        flatTopThreshold = exp(-3)*exp(-3);
        coneApertureProfile(coneApertureProfile>=flatTopThreshold) = flatTopThreshold;
        coneApertureProfile(coneApertureProfile<flatTopThreshold) = 0;
        coneApertureProfile = coneApertureProfile / flatTopThreshold;
    else  
        coneApertureProfile = exp(-(X/coneCharacteristicRadiusDegs).^2) .* exp(-(Y/coneCharacteristicRadiusDegs).^2);
    end
    
    % cone aperture sampling grid
    deltaFunctionMap2D = zeros(length(thePSFsupportDegs), length(thePSFsupportDegs));
    for iCone = 1:size(conePosDegs,1)
        coneXpos = conePosDegs(iCone,1);
        coneYpos = conePosDegs(iCone,2);
        [~,ix] = min(abs(coneXpos - thePSFsupportDegs));
        [~,iy] = min(abs(coneYpos - thePSFsupportDegs));
        if ( (coneXpos>=thePSFsupportDegs(1)) && (coneXpos<=thePSFsupportDegs(end)) && ...
             (coneYpos>=thePSFsupportDegs(1)) && (coneYpos<=thePSFsupportDegs(end)) )
            deltaFunctionMap2D(iy,ix) = 1;
        end
    end
    
    % Cone aperture sampling map
    coneSamplingMap = conv2(deltaFunctionMap2D, coneApertureProfile, 'same');
    
    % Fine the 6 cones surrounding the central-most cone
    % (7 total = 6 neighboring + central cone)
    [~, idx] = sort(sqrt(sum(conePosDegs.^2,2)), 'ascend');
    neigboringConeIndices = idx(1:7);
    
    % Generate gaussian pooling configs
    gaussianPoolingConfigs = gaussianPoolingCenterForNconeInputRGC(coneCharacteristicRadiusDegs, conePosDegs);
    
    % Deconvolution data is a container
    deconvolutionStruct.data = containers.Map();  
    
    % Whether to visualize all the intermediate maps
    visualizeIntermediateMaps = true;
    
    % Analyze each pooling config
    for poolingConfigIndex = 1:numel(gaussianPoolingConfigs)
        
        fprintf('Analyzing %d-cone input center subregion.\n',  gaussianPoolingConfigs{poolingConfigIndex}.coneInputsNum);
        inputConeIndicesForAllCombinations = gaussianPoolingConfigs{poolingConfigIndex}.inputConeIndices;
        
        % Preallocate memory
        conesNum = size(conePosDegs,1);
        continuousPoolingMapsInRetinalSpace = zeros(size(inputConeIndicesForAllCombinations,1), size(X,1), size(X,2));
        continuousPoolingMapsInVisualSpace = zeros(size(inputConeIndicesForAllCombinations,1), size(X,1), size(X,2));
        withinConePoolingMapsInRetinalSpace = zeros(size(inputConeIndicesForAllCombinations,1), size(X,1), size(X,2));
        withinConePoolingMapsInVisualSpace = zeros(size(inputConeIndicesForAllCombinations,1), size(X,1), size(X,2));
        integratedConePoolingMapsInRetinalSpace = zeros(size(inputConeIndicesForAllCombinations,1),conesNum);
        integratedConePoolingMapsInVisualSpace = zeros(size(inputConeIndicesForAllCombinations,1),conesNum);
    
        neigboringConeWeightsInRetinalSpace = zeros(size(inputConeIndicesForAllCombinations,1), numel(neigboringConeIndices));
        neigboringConeWeightsInVisualSpace = zeros(size(inputConeIndicesForAllCombinations,1), numel(neigboringConeIndices));
        
        for inputConeCombination = 1:size(inputConeIndicesForAllCombinations,1)
            inputConeIndices = inputConeIndicesForAllCombinations(inputConeCombination,:);
        
            deltaFunctionMap2D = zeros(length(thePSFsupportDegs), length(thePSFsupportDegs));
            for iCone = 1:numel(inputConeIndices)
                coneXpos = conePosDegs(inputConeIndices(iCone),1);
                coneYpos = conePosDegs(inputConeIndices(iCone),2);
                [~,ix] = min(abs(coneXpos - thePSFsupportDegs));
                [~,iy] = min(abs(coneYpos - thePSFsupportDegs));
                if ( (coneXpos>=thePSFsupportDegs(1)) && (coneXpos<=thePSFsupportDegs(end)) && ...
                     (coneYpos>=thePSFsupportDegs(1)) && (coneYpos<=thePSFsupportDegs(end)) )
                    deltaFunctionMap2D(iy,ix) = 1;
                end
            end

            % Cone aperture sampling map
            continuousPoolingMapInRetinalSpace = conv2(deltaFunctionMap2D, coneApertureProfile, 'same');
                      
            % The continous cone pooling in visual space
            continuousPoolingMapInVisualSpace = conv2(continuousPoolingMapInRetinalSpace, thePSF, 'same');
        
            % Pooling weights within cone apertures mapped in retinal space
            withinConePoolingMapInRetinalSpace = coneSamplingMap .* continuousPoolingMapInRetinalSpace;
            withinConePoolingMapInRetinalSpace(withinConePoolingMapInRetinalSpace>0.01*max(withinConePoolingMapInRetinalSpace(:))) = 1;
            
            % Pooling weights within cone apertures mapped in visual space
            withinConePoolingMapInVisualSpace =  coneSamplingMap .* continuousPoolingMapInVisualSpace;
            
            % Integrated cone pooling weights
            for iCone = 1:conesNum
                xr = conePosDegs(iCone,1) + coneAperturesDegs(iCone)*0.5*cosd(0:20:360);
                yr = conePosDegs(iCone,2) + coneAperturesDegs(iCone)*0.5*sind(0:20:360);
                [in,on] = inpolygon(X(:), Y(:), xr, yr);
                allIndices = find((in==true));
                if (numel(allIndices) == 0)
                    integratedConePoolingMapsInRetinalSpace(inputConeCombination, iCone) = 0;
                    integratedConePoolingMapsInVisualSpace(inputConeCombination, iCone) = 0;
                else
                    integratedConePoolingMapsInRetinalSpace(inputConeCombination, iCone) = mean(withinConePoolingMapInRetinalSpace(allIndices));
                    integratedConePoolingMapsInVisualSpace(inputConeCombination, iCone)  = mean(withinConePoolingMapInVisualSpace(allIndices));
                    if (integratedConePoolingMapsInVisualSpace(inputConeCombination, iCone)<0) || (integratedConePoolingMapsInVisualSpace(inputConeCombination, iCone)>1)
                        fprintf('How can this be. Integreated cone pooling is R%f\n', integratedConePoolingMapsInVisualSpace(inputConeCombination, iCone))
                    end
                end
            end
        
            continuousPoolingMapsInRetinalSpace(inputConeCombination,:,:) = continuousPoolingMapInRetinalSpace;
            continuousPoolingMapsInVisualSpace(inputConeCombination,:,:) = continuousPoolingMapInVisualSpace;
            withinConePoolingMapsInRetinalSpace(inputConeCombination,:,:) = withinConePoolingMapInRetinalSpace;
            withinConePoolingMapsInVisualSpace(inputConeCombination,:,:) = withinConePoolingMapInVisualSpace;
        
            neigboringConeWeightsInRetinalSpace(inputConeCombination,:) = integratedConePoolingMapsInRetinalSpace(inputConeCombination,neigboringConeIndices);
            neigboringConeWeightsInVisualSpace(inputConeCombination,:) = integratedConePoolingMapsInVisualSpace(inputConeCombination,neigboringConeIndices);
        end %  inputConeCombination
        

        % Fit integratedConePoolingMapInRetinalSpace to extract retinal pooling sigmas and gains
        [fittedRetinalActivationMaps, fittedRetinalSlices, rfSigmasRetinal, rfGainRetinal, hiResPSFsupportDegs] = ...
             fitConePoolingMaps(integratedConePoolingMapsInRetinalSpace, ...
             coneCharacteristicRadiusDegs, inputConeIndicesForAllCombinations, ...
             conePosDegs(:,1), conePosDegs(:,2), thePSFsupportDegs);

        % Fit integratedConePoolingMapInVisualSpace to extract visual pooling sigmas and gains
        [fittedVisualActivationMaps, fittedVisualSlices, rfSigmasVisual, rfGainVisual] = ...
             fitConePoolingMaps(integratedConePoolingMapsInVisualSpace, ...
             coneCharacteristicRadiusDegs, [], ...
             conePosDegs(:,1), conePosDegs(:,2), thePSFsupportDegs);

        % Extracted data
        keyName = sprintf('%d-coneInput',gaussianPoolingConfigs{poolingConfigIndex}.coneInputsNum);
        deconvolutionStruct.data(keyName) = struct(...
            'visualGainAttenuation', min(rfGainVisual./rfGainRetinal), ...
            'minRetinalSigma', mean(min(rfSigmasRetinal,[],2)), ...
            'maxRetinalSigma', mean(max(rfSigmasRetinal,[],2)), ...
            'minVisualSigma', mean(min(rfSigmasVisual,[],2)), ...
            'maxVisualSigma', mean(max(rfSigmasVisual,[],2)) ...
        );
        
        if (visualizeIntermediateMaps)
            % Put all intermediate maps in a struct
            allMaps = struct;
            allMaps.rfSigmasRetinal = rfSigmasRetinal;
            allMaps.rfGainRetinal = rfGainRetinal;
            allMaps.rfSigmasVisual = rfSigmasVisual;
            allMaps.rfGainVisual = rfGainVisual;
            allMaps.integratedConePoolingMapsInRetinalSpace = integratedConePoolingMapsInRetinalSpace;
            allMaps.integratedConePoolingMapsInVisualSpace = integratedConePoolingMapsInVisualSpace;
            allMaps.withinConePoolingMapsInVisualSpace = withinConePoolingMapsInVisualSpace;
            allMaps.withinConePoolingMapsInRetinalSpace = withinConePoolingMapsInRetinalSpace;
            allMaps.continuousPoolingMapsInVisualSpace = continuousPoolingMapsInVisualSpace;
            allMaps.continuousPoolingMapsInRetinalSpace = continuousPoolingMapsInRetinalSpace;
            allMaps.neigboringConeWeightsInVisualSpace = neigboringConeWeightsInVisualSpace;
            allMaps.neigboringConeWeightsInRetinalSpace = neigboringConeWeightsInRetinalSpace;
            allMaps.fittedRetinalActivationMaps = fittedRetinalActivationMaps;
            allMaps.fittedVisualActivationMaps = fittedVisualActivationMaps;
            allMaps.fittedRetinalSlices = fittedRetinalSlices;
            allMaps.fittedVisualSlices = fittedVisualSlices;
            % Visualize everything
            visualizeIntemediateDeconvolutionAnalysisMaps(PolansSubjectID, patchEccDegs, ...
                gaussianPoolingConfigs{poolingConfigIndex}.coneInputsNum, ...
                conePosDegs, coneAperturesDegs, thePSF, thePSFsupportDegs, hiResPSFsupportDegs, ...
                allMaps, coneCharacteristicRadiusDegs);
        end % visualizeIntermediateMaps
    end % poolingConfigIndex   
end
    

function visualizeIntemediateDeconvolutionAnalysisMaps(PolansSubjectID, patchEccDegs, ...
    coneInputsNum, conePosDegs, coneAperturesDegs, thePSF, thePSFsupportDegs, hiResPSFsupportDegs, ...
    allMaps, coneCharacteristicRadiusDegs)
    
    visualizeSummaryMaps(100+coneInputsNum, coneInputsNum, ...
        allMaps.neigboringConeWeightsInRetinalSpace, allMaps.neigboringConeWeightsInVisualSpace, ...
        allMaps.rfSigmasVisual, allMaps.rfSigmasRetinal, allMaps.rfGainVisual, allMaps.rfGainRetinal, ...
         coneCharacteristicRadiusDegs);
    
    if (1==2)
        visualizeDetailedMaps(1000*coneInputsNum, coneInputsNum, ...
            PolansSubjectID, patchEccDegs, ...
            conePosDegs, coneAperturesDegs, thePSF, thePSFsupportDegs, hiResPSFsupportDegs, ...
            allMaps.continuousPoolingMapsInRetinalSpace, allMaps.continuousPoolingMapsInVisualSpace, ...
            allMaps.withinConePoolingMapsInRetinalSpace, allMaps.withinConePoolingMapsInVisualSpace, ...
            allMaps.integratedConePoolingMapsInRetinalSpace, allMaps.integratedConePoolingMapsInVisualSpace, ...
            allMaps.fittedRetinalActivationMaps, allMaps.fittedVisualActivationMaps, ...
            allMaps.fittedRetinalSlices, allMaps.fittedVisualSlices, ...
            coneCharacteristicRadiusDegs);
    end
end

function visualizeSummaryMaps(figNo, coneInputsNum, ...
    neigboringConeWeightsInRetinalSpace, neigboringConeWeightsInVisualSpace, ...
    rfSigmasVisual, rfSigmasRetinal, rfGainVisual, rfGainRetinal, ...
    coneCharacteristicRadiusDegs)

    % Summary figure showing the amplitudes and sigmas for the
    % different Gaussian locations
    hFig = figure(figNo); clf;
    set(hFig, 'Position', [100 100 26 13], 'Name', sprintf('Summary figure for %d-cone input scenario', coneInputsNum));

    theAxesGrid = plotlab.axesGrid(hFig, ...
            'leftMargin', 0.04, ...
            'bottomMargin', 0.04, ...
            'rightMargin', 0.01, ...
            'widthMargin', 0.05, ...
            'heightMargin', 0.08, ...
            'topMargin', 0.03, ...
            'rowsNum', 2, ...
            'colsNum', 5);

    % Plot the cone weights for the neighboring cones
    for iAngle = 1:size(neigboringConeWeightsInRetinalSpace,1)
        switch(iAngle)
            case 1
                row = 1; col = 3;
            case 2
                row = 1; col = 4;
            case 3
                row = 1; col = 5;
            case 4
                row = 2; col = 3;
            case 5
                row = 2; col = 4;
            case 6
                row = 2; col = 5;    
        end

        ax = theAxesGrid{row,col};
        weights = [squeeze(neigboringConeWeightsInRetinalSpace(iAngle,:)); squeeze(neigboringConeWeightsInVisualSpace(iAngle,:))];
        bar(ax, 1:size(neigboringConeWeightsInRetinalSpace,2), weights, 1.0);
        set(ax, 'XLim', [0 8], 'XTick', 1:7, 'XTickLabel', {'c', 's1', 's2', 's3', 's4', 's5', 's6'}, 'YLim', [0 1]);
        xlabel('cone id');
        ylabel('pooling weight');
    end

    % Sigma range
    maxSigma = max([max(rfSigmasVisual(:)) max(rfSigmasRetinal(:))]);
    minSigma = min([min(rfSigmasVisual(:)) min(rfSigmasRetinal(:))]);
    dSigma = (maxSigma-minSigma)/4;
    sigmaRange = ([minSigma maxSigma] + dSigma*[-1 1]);

    % Fitted sigmas (min,max) of the best 2D elliptical Gaussian in retinal space
    ax = theAxesGrid{1,1}; 
    visualizeFittedGaussianSigmas(ax, rfSigmasRetinal, sigmaRange,  coneCharacteristicRadiusDegs, 'retinal');

    % Fitted sigmas (min,max) of the best 2D elliptical Gaussian in visual space
    ax = theAxesGrid{2,1}; 
    visualizeFittedGaussianSigmas(ax, rfSigmasVisual, sigmaRange,  coneCharacteristicRadiusDegs, 'visual');

     % Fitted gain of the best 2D elliptical Gaussian in retinal space
    ax = theAxesGrid{1,2};
    visualizeFittedGaussianGain(ax,rfGainRetinal);
    
    % Fitted gain of the best 2D elliptical Gaussian in visual space
    ax = theAxesGrid{2,2};
    visualizeFittedGaussianGain(ax,rfGainVisual);
end
    
function visualizeDetailedMaps(baseFigNo, coneInputsNum, PolansSubjectID, patchEccDegs, ...
    conePosDegs, coneAperturesDegs, thePSF, thePSFsupportDegs, hiResPSFsupportDegs, ...
    continuousPoolingMapsInRetinalSpace, continuousPoolingMapsInVisualSpace, ...
    withinConePoolingMapsInRetinalSpace, withinConePoolingMapsInVisualSpace, ...
    integratedConePoolingMapsInRetinalSpace, integratedConePoolingMapsInVisualSpace, ...
    fittedRetinalActivationMaps, fittedVisualActivationMaps, ...
    fittedRetinalSlices, fittedVisualSlices, ...
    coneCharacteristicRadiusDegs)

    % Visualized cone range
    visualizecConesNumAcross = 6;
    visualizeSpaceRange = visualizecConesNumAcross * coneCharacteristicRadiusDegs*3*2*[-1 1];
    visualizedSpaceTicksDegs = -0.2:0.1:0.2;

    % Plot detail analyses for all Gaussian locations
    for iAngle = 1:size(continuousPoolingMapsInRetinalSpace,1)

        continuousPoolingMapInRetinalSpace = squeeze(continuousPoolingMapsInRetinalSpace(iAngle,:,:));
        continuousPoolingMapInVisualSpace = squeeze(continuousPoolingMapsInVisualSpace(iAngle,:,:));

        withinConePoolingMapInRetinalSpace = squeeze(withinConePoolingMapsInRetinalSpace(iAngle,:,:));
        withinConePoolingMapInVisualSpace = squeeze(withinConePoolingMapsInVisualSpace(iAngle,:,:));

        integratedConePoolingMapInRetinalSpace = squeeze(integratedConePoolingMapsInRetinalSpace(iAngle,:));
        integratedConePoolingMapInVisualSpace = squeeze(integratedConePoolingMapsInVisualSpace(iAngle,:));

        fittedRetinalActivationMap = squeeze(fittedRetinalActivationMaps(iAngle,:,:));
        fittedVisualActivationMap = squeeze(fittedVisualActivationMaps(iAngle,:,:));

        fittedRetinalSlice = squeeze(fittedRetinalSlices(iAngle,:));
        fittedVisualSlice = squeeze(fittedVisualSlices(iAngle,:));

        hFig = figure(baseFigNo+iAngle); clf;
        set(hFig, 'Position', [100 100 26 13], 'Name', sprintf('Detailed figure for %d-cone input scenario, angle index %d', coneInputsNum, iAngle));

        theAxesGrid = plotlab.axesGrid(hFig, ...
            'leftMargin', 0.04, ...
            'bottomMargin', 0.04, ...
            'rightMargin', 0.01, ...
            'widthMargin', 0.05, ...
            'heightMargin', 0.08, ...
            'topMargin', 0.03, ...
            'rowsNum', 2, ...
            'colsNum', 4);

        % Plot the PSF and the cone apertures
        ax = theAxesGrid{2,1};
        coneIndicesToPlot = find(...
            abs(conePosDegs(:,1)) <= max(visualizeSpaceRange) & ...
            abs(conePosDegs(:,2)) <= max(visualizeSpaceRange) );

        visualineConeAperturesAndPSF(ax, conePosDegs(coneIndicesToPlot,:), coneAperturesDegs(coneIndicesToPlot,:), ...
            thePSF, thePSFsupportDegs, visualizeSpaceRange, ...
            'xTicks', [], 'yTicks', visualizedSpaceTicksDegs, ...
            'titleLabel', sprintf('cone mosaic and PSF\n Subject %d @(%2.1f, %2.1f) degrees', PolansSubjectID, patchEccDegs(1), patchEccDegs(2)));

        % Plot the continuous retinal cone pooling weights
        ax = theAxesGrid{1,2};
        visualizeContinuousPoolingAndPSF(ax, thePSFsupportDegs, visualizeSpaceRange,...
            conePosDegs(coneIndicesToPlot,:), coneAperturesDegs(coneIndicesToPlot,:),  ...
            continuousPoolingMapInRetinalSpace, [],  ...
            'xTicks', [], 'yTicks', visualizedSpaceTicksDegs, ...
            'titleLabel', 'retinally-mapped cone pooling');

        % Plot the continuous visual cone pooling weights and the PSF
        ax = theAxesGrid{2,2};
        visualizeContinuousPoolingAndPSF(ax, thePSFsupportDegs, visualizeSpaceRange, ...
            conePosDegs(coneIndicesToPlot,:), coneAperturesDegs(coneIndicesToPlot,:), ...
            continuousPoolingMapInVisualSpace, thePSF, ...
            'xTicks', visualizedSpaceTicksDegs, 'yTicks', visualizedSpaceTicksDegs, ...
            'titleLabel', 'visually-mapped pooling');
        
        % Plot the within-aperture retinal cone pooling weights
        ax = theAxesGrid{1,3};
        visualizeWithinConeAperturePooling(ax, thePSFsupportDegs, visualizeSpaceRange, ...
            conePosDegs(coneIndicesToPlot,:), coneAperturesDegs(coneIndicesToPlot,:), ...
            withinConePoolingMapInRetinalSpace, ...
            'xTicks', [], 'yTicks', [], ...
            'titleLabel','retinally-mapped pooling');

        % Plot the within-aperture visual cone pooling weights
        ax = theAxesGrid{2,3};
        visualizeWithinConeAperturePooling(ax, thePSFsupportDegs, visualizeSpaceRange, ...
            conePosDegs(coneIndicesToPlot,:), coneAperturesDegs(coneIndicesToPlot,:), ...
            withinConePoolingMapInVisualSpace, ...
            'xTicks', visualizedSpaceTicksDegs, 'yTicks', [], ...
            'titleLabel', 'visually-mapped pooling');

        % Plot the integrated retinal cone pooling weights, along with the best-fit 2D elliptical Gaussian profile
        ax = theAxesGrid{1,4};
        visualizeIntegratedConePooling(ax, thePSFsupportDegs, visualizeSpaceRange, ...
            conePosDegs(coneIndicesToPlot,:), coneAperturesDegs(coneIndicesToPlot,:), ...
            integratedConePoolingMapInRetinalSpace(coneIndicesToPlot), fittedRetinalActivationMap, hiResPSFsupportDegs, ...
            'xTicks', [], 'yTicks', [], ...
            'titleLabel','retinally-mapped pooling');

        % Plot the integrated visual cone pooling weights, along with the best-fit 2D elliptical Gaussian profile
        ax = theAxesGrid{2,4};
        visualizeIntegratedConePooling(ax, thePSFsupportDegs, visualizeSpaceRange, ...
            conePosDegs(coneIndicesToPlot,:), coneAperturesDegs(coneIndicesToPlot,:), ...
            integratedConePoolingMapInVisualSpace(coneIndicesToPlot),  fittedVisualActivationMap, hiResPSFsupportDegs, ...
            'xTicks', visualizedSpaceTicksDegs, 'yTicks', [], ...
            'titleLabel','visually-mapped pooling');

        % Plot the activation of cones (in retinal and visual space view) along the horizontal meridian
        ax = theAxesGrid{1,1};
        midRowCones = find(abs(conePosDegs(:,2)) < 1e-3);
        hold(ax,'on')
        for iCone = 1:numel(midRowCones)
           scatter(ax, conePosDegs(midRowCones,1), integratedConePoolingMapInRetinalSpace(midRowCones), 'o', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', [0 0 0]);
           scatter(ax, conePosDegs(midRowCones,1), integratedConePoolingMapInVisualSpace(midRowCones),  'o', 'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', [1 0 0]);
        end
        plot(ax, hiResPSFsupportDegs, fittedRetinalSlice, 'k-');
        plot(ax, hiResPSFsupportDegs, fittedVisualSlice, 'r-');
        axis(ax, 'square');
        set(ax, 'XLim', visualizeSpaceRange, 'YLim', [0 1]);
        colormap(brewermap(1024, '*YlGnBu'));
        drawnow;
    end % iAngle
       
 end 



function gaussianPoolingCenters = gaussianPoolingCenterForNconeInputRGC(coneCharacteristicRadiusDegs, conePosDegs)

    % For 1-cone input RGCs, place the center of the Gaussian at [0 0]
    gaussianPoolingCenters{1} = struct;
    gaussianPoolingCenters{1}.coneInputsNum = 1; 
    gaussianPoolingCenters{1}.poolingCenters = [0 0];
    [~,inputConeIndicesForThisAngle] = pdist2(conePosDegs, gaussianPoolingCenters{1}.poolingCenters, 'euclidean', 'smallest', 1);
    gaussianPoolingCenters{1}.inputConeIndices = inputConeIndicesForThisAngle;
    
    
    % For 2-cone input RGCs, generate 6 Gaussians, exactly midway between 
    % the center cone and its neigboring neighbor
    poolingCenters = zeros(6,1);
    inputConeIndices = zeros(6,2);
    for iAngle = 1:size(poolingCenters,1)
        offsetRadius = coneCharacteristicRadiusDegs/WatsonRGCModel.coneApertureToDiameterRatio*3;
        poolingCenters(iAngle,1) = offsetRadius*cosd(30+iAngle*60);
        poolingCenters(iAngle,2) = offsetRadius*sind(30+iAngle*60);
        [~,inputConeIndicesForThisAngle] = pdist2(conePosDegs, poolingCenters(iAngle,:), 'euclidean', 'smallest', 2);
        inputConeIndices(iAngle,:) = inputConeIndicesForThisAngle;
    end
    gaussianPoolingCenters{2} = struct;
    gaussianPoolingCenters{2}.coneInputsNum = 2; 
    gaussianPoolingCenters{2}.poolingCenters = poolingCenters;
    gaussianPoolingCenters{2}.inputConeIndices = inputConeIndices;
    
    % For 3-cone input RGCs, generate 6 Gaussians, exactly midway between 
    % the center cone and its neigboring 2 neighbors
    poolingCenters = zeros(6,1);
    inputConeIndices = zeros(6,3);
    for iAngle = 1:size(poolingCenters,1)
        offsetRadius = 2/sqrt(3.0)*coneCharacteristicRadiusDegs/WatsonRGCModel.coneApertureToDiameterRatio*3;
        poolingCenters(iAngle,1) = offsetRadius*cosd(iAngle*60);
        poolingCenters(iAngle,2) = offsetRadius*sind(iAngle*60);
        [~,inputConeIndicesForThisAngle] = pdist2(conePosDegs, poolingCenters(iAngle,:), 'euclidean', 'smallest', 3);
        inputConeIndices(iAngle,:) = inputConeIndicesForThisAngle;
    end
    gaussianPoolingCenters{3} = struct;
    gaussianPoolingCenters{3}.coneInputsNum = 3; 
    gaussianPoolingCenters{3}.poolingCenters = poolingCenters;
    gaussianPoolingCenters{3}.inputConeIndices = inputConeIndices;
    
    % For 4-cone input RGCs, generate 8 Gaussians,
    poolingCenters = zeros(8,1);
    inputConeIndices = zeros(8,4);
    for iAngle = 1:size(poolingCenters,1)
        offsetRadius = 2/sqrt(3.0)*coneCharacteristicRadiusDegs/WatsonRGCModel.coneApertureToDiameterRatio*3;
        poolingCenters(iAngle,1) = offsetRadius*cosd((iAngle-1)*45);
        poolingCenters(iAngle,2) = offsetRadius*sind((iAngle-1)*45);
        [~,inputConeIndicesForThisAngle] = pdist2(conePosDegs, poolingCenters(iAngle,:), 'euclidean', 'smallest', 4);
        inputConeIndices(iAngle,:) = inputConeIndicesForThisAngle;
    end
    gaussianPoolingCenters{4} = struct;
    gaussianPoolingCenters{4}.coneInputsNum = 4; 
    gaussianPoolingCenters{4}.poolingCenters = poolingCenters;
    gaussianPoolingCenters{4}.inputConeIndices = inputConeIndices;
    
    % For 5 and 6cone input RGCs, generate 2 Gaussians
    for coneInputsNum = 5:5
        poolingCenters = zeros(2,1);
        inputConeIndices = zeros(2,coneInputsNum);
        for iAngle = 1:size(poolingCenters,1)
            offsetRadius = coneCharacteristicRadiusDegs/WatsonRGCModel.coneApertureToDiameterRatio*3;
            poolingCenters(iAngle,1) = offsetRadius*cosd((iAngle-1)*180);
            poolingCenters(iAngle,2) = offsetRadius*sind((iAngle-1)*180);
            [~,inputConeIndicesForThisAngle] = pdist2(conePosDegs, poolingCenters(iAngle,:), 'euclidean', 'smallest', coneInputsNum);
            inputConeIndices(iAngle,:) = inputConeIndicesForThisAngle;
        end
        gaussianPoolingCenters{coneInputsNum} = struct;
        gaussianPoolingCenters{coneInputsNum}.coneInputsNum = coneInputsNum; 
        gaussianPoolingCenters{coneInputsNum}.poolingCenters = poolingCenters;
        gaussianPoolingCenters{coneInputsNum}.inputConeIndices = inputConeIndices;
    end
    
    % For 7+cone input RGCs, generate 6 Gaussians, exactly midway between 
    % the center cone and its neigboring 2 neighbors
    for coneInputsNum = 6:100
        gaussianPoolingCenters{coneInputsNum} = struct;
        gaussianPoolingCenters{coneInputsNum}.coneInputsNum = coneInputsNum; 
        gaussianPoolingCenters{coneInputsNum}.poolingCenters = [0 0];
        [~,inputConeIndices] = pdist2(conePosDegs, [0 0], 'euclidean', 'smallest', coneInputsNum);
        gaussianPoolingCenters{coneInputsNum}.inputConeIndices(1,:) = inputConeIndices;
    end
    
end

function [fittedActivationMapsFull2D, fittedActivationSlices,  rfSigmas, rfGain, hiResPSFsupportDegs] = fitConePoolingMaps(...
    conePoolingMaps, minCharacteristicRadiusDegs, inputConeIndicesForAllCombinations, X, Y, thePSFsupportDegs)


    for inputConfigIndex = 1:size(conePoolingMaps,1)
        if (~isempty(inputConeIndicesForAllCombinations))
            inputConeIndices = inputConeIndicesForAllCombinations(inputConfigIndex,:)
        else
            inputConeIndices = [];
        end
        [tmp1, tmp2, tmp3, tmp4, hiResPSFsupportDegs] = fitActivationMap(...
            squeeze(conePoolingMaps(inputConfigIndex,:)), minCharacteristicRadiusDegs, ...
            inputConeIndices, X, Y, thePSFsupportDegs);
        fittedActivationMapsFull2D(inputConfigIndex,:,:) = tmp1;
        fittedActivationSlices(inputConfigIndex,:) = tmp2;  ...
        rfSigmas(inputConfigIndex,:) = tmp3;
        rfGain(inputConfigIndex) = tmp4;
    end
    
end

function [fittedActivationMapFull2D, fittedActivationSlice,  rfSigmas, rfGain, hiResPSFsupportDegs] = fitActivationMap(...
    activationMap, minCharacteristicRadiusDegs, inputConeIndices, X, Y, thePSFsupportDegs)

    deltaX = thePSFsupportDegs(2)-thePSFsupportDegs(1);
    hiResPSFsupportDegs = min(thePSFsupportDegs):deltaX/4:max(thePSFsupportDegs);
    
    [fittedParams, rfFunction] = ...
        fitElliptical2DGausianToRF(X(:), Y(:), activationMap(:), deltaX/10, minCharacteristicRadiusDegs, inputConeIndices, [0 0]);
    
    rfGain = fittedParams(1);
    if (numel(fittedParams) == 4)
        rfSigmas = fittedParams(4)*[1 1];
    else
        rfSigmas = [fittedParams(4) fittedParams(5)];
    end
    
    [XXX,YYY] = meshgrid(hiResPSFsupportDegs, hiResPSFsupportDegs);
    xyData(:,1) = XXX(:);
    xyData(:,2) = YYY(:);
    
    fittedActivationMap = rfFunction(fittedParams,xyData);
    fittedActivationMapFull2D = reshape(fittedActivationMap, [numel(hiResPSFsupportDegs) numel(hiResPSFsupportDegs)]);
    
    [~,midRow] = min(abs(hiResPSFsupportDegs));
    fittedActivationSlice = squeeze(fittedActivationMapFull2D(midRow,:));
end


function visualizeFittedGaussianSigmas(ax, rfSigmas, sigmaRange,  coneCharacteristicRadiusDegs, sigmaDomain)
    configIDs = 1:size(rfSigmas,1);
    scatter(ax,configIDs, min(rfSigmas,[],2)*60, 'o', 'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', [1 0.3 0.3]); 
    hold(ax, 'on');
    scatter(ax,configIDs, max(rfSigmas,[],2)*60, 'o', 'MarkerFaceColor', [0.5 0.3 0.3], 'MarkerEdgeColor', [0.5 0. 0.]); 
    line(ax, [0 8], coneCharacteristicRadiusDegs*60*[1 1], 'LineStyle','--', 'Color', [0 0 0], 'LineWidth', 1.5);
    sigmaRange = sigmaRange*60;
    if (sigmaRange(1)>10)
        sigmaRange = round(sigmaRange);
    end
    set(ax, 'YLim', sigmaRange);
    set(ax, 'XLim', [0 8], 'XTick', 1:7, 'XTickLabel', {'c', 's1', 's2', 's3', 's4', 's5', 's6'}, 'YLim', [0 1]);
    xlabel('configuration');
    ylabel(ax, sprintf('%s characteristic radius (arc min)', sigmaDomain));
    box(ax, 'on'); 
    
end

function visualizeFittedGaussianGain(ax,rfGainVisual)
    configIDs = 1:numel(rfGainVisual);
    scatter(ax, configIDs, rfGainVisual, 'd', 'MarkerFaceColor', [0.5 0.5 1], 'MarkerEdgeColor', [0 0 1]);
    ylabel(ax, 'visual gain');
    set(ax, 'XLim', [0 8], 'XTick', 1:7, 'XTickLabel', {'c', 's1', 's2', 's3', 's4', 's5', 's6'}, 'YLim', [0 1]);
    xlabel('configuration');
    box(ax, 'on');    
end

        
function visualineConeAperturesAndPSF(ax, conePosDegs, coneAperturesDegs, thePSF, thePSFsupportDegs, visualizeSpaceRange, varargin)
    % Parse input
    p = inputParser;
    p.addParameter('xTicks', [], @isnumeric);
    p.addParameter('yTicks', [], @isnumeric);
    p.addParameter('titleLabel', '', @ischar);
    p.parse(varargin{:});
   
    cmap = parula(1024);
    [X,Y] = meshgrid(thePSFsupportDegs,thePSFsupportDegs);
    edgeColor = [0.5 0.5 0.5]; faceAlpha = 1.0;
    imagesc(ax, thePSFsupportDegs, thePSFsupportDegs, zeros(numel(thePSFsupportDegs), numel(thePSFsupportDegs)));
    hold(ax, 'on')
    plotConeApertures(ax, conePosDegs, coneAperturesDegs, zeros(size(coneAperturesDegs)), edgeColor, faceAlpha); hold on
    contour(ax, X,Y, thePSF/max(thePSF(:)), 0:0.2:1, 'LineColor', 'c', 'LineWidth', 1.5);
    colormap(ax, cmap);
    
    axis(ax, 'xy'); axis(ax, 'square'); box(ax, 'on');
    set(ax, 'XLim', visualizeSpaceRange, 'YLim', visualizeSpaceRange, 'CLim', [0 1]);
    set(ax, 'XTick', p.Results.xTicks, 'YTick',  p.Results.xTicks);
    if (~isempty(p.Results.titleLabel))
        title(ax, p.Results.titleLabel);
    end
    xlabel(ax, 'space (degs)');

end

function visualizeContinuousPoolingAndPSF(ax, thePSFsupportDegs, visualizeSpaceRange, conePosDegs, coneAperturesDegs, continuousActivation, thePSF, ...
        varargin)
    
    % Parse input
    p = inputParser;
    p.addParameter('xTicks', [], @isnumeric);
    p.addParameter('yTicks', [], @isnumeric);
    p.addParameter('titleLabel', '', @ischar);
    p.parse(varargin{:});
    
    [~,idx] = max(continuousActivation(:));
    [midRow, midCol] = ind2sub(size(continuousActivation),idx);
    activationSlice = squeeze(continuousActivation(midRow,:));
    offset = visualizeSpaceRange(1);
    gain =  (visualizeSpaceRange(2)-visualizeSpaceRange(1));
    
    hold(ax, 'on');
    imagesc(ax, thePSFsupportDegs, thePSFsupportDegs, continuousActivation, [0 1]);
    edgeColor = [0.5 0.5 0.5]; faceAlpha = 0.1;
    plotConeApertures(ax, conePosDegs, coneAperturesDegs, zeros(size(coneAperturesDegs)), edgeColor, faceAlpha);
    line(ax, thePSFsupportDegs, offset+gain*activationSlice, 'Color', [1 1 0], 'LineWidth', 1.5);
    
%     if (~isempty(thePSF))
%         [X,Y] = meshgrid(thePSFsupportDegs, thePSFsupportDegs);
%         contour(ax, X,Y, thePSF/max(thePSF(:)), 0:0.2:1, 'LineColor', 'r', 'LineWidth', 1.5);
%     end
    axis(ax, 'xy'); axis(ax,'square'); box(ax, 'on');
    set(ax, 'XLim', visualizeSpaceRange, 'YLim', visualizeSpaceRange, 'CLim', [0 1]);
    
    set(ax, 'XTick', p.Results.xTicks, 'YTick', p.Results.yTicks);
    if (~isempty(p.Results.titleLabel))
        title(ax, sprintf('%s\n(continuous)', p.Results.titleLabel));
    end

end

function visualizeWithinConeAperturePooling(ax, thePSFsupportDegs, visualizeSpaceRange, conePosDegs, coneAperturesDegs, coneBasedSubgegionActivation, varargin)
    % Parse input
    p = inputParser;
    p.addParameter('xTicks', [], @isnumeric);
    p.addParameter('yTicks', [], @isnumeric);
    p.addParameter('titleLabel', '', @ischar);
    p.parse(varargin{:});
    
%     [~,idx] = max(coneBasedSubgegionActivation(:));
%     [midRow, midCol] = ind2sub(size(coneBasedSubgegionActivation),idx);
%     activationSlice = squeeze(coneBasedSubgegionActivation(midRow,:));
%     offset = visualizeSpaceRange(1);
%     gain =  (visualizeSpaceRange(2)-visualizeSpaceRange(1));
    
    imagesc(ax, thePSFsupportDegs, thePSFsupportDegs, coneBasedSubgegionActivation, [0 1]); hold on;
    edgeColor = [0.5 0.5 0.5]; faceAlpha = 0.1;
    plotConeApertures(ax, conePosDegs, coneAperturesDegs, zeros(size(coneAperturesDegs)), edgeColor, faceAlpha);
    %line(ax, thePSFsupportDegs, offset+gain*activationSlice, 'Color', [1 1 0], 'LineWidth', 1.5);
    
    set(ax, 'XLim', visualizeSpaceRange, 'YLim', visualizeSpaceRange);
    axis(ax, 'xy');
    axis(ax, 'square');
    set(ax, 'XLim', visualizeSpaceRange, 'YLim', visualizeSpaceRange, 'CLim', [0 1]);
    
    set(ax, 'XTick', p.Results.xTicks, 'YTick', p.Results.yTicks); box(ax, 'on');
    if (~isempty(p.Results.titleLabel))
        title(ax, sprintf('%s\n(within aperture)', p.Results.titleLabel));
    end

end

function visualizeIntegratedConePooling(ax, thePSFsupportDegs, visualizeSpaceRange, conePosDegs, coneAperturesDegs, ...
    integratedConeBasedSubgegionActivation, fittedActivationMap, hiResPSFsupportDegs, varargin)
    % Parse input
    p = inputParser;
    p.addParameter('xTicks', [], @isnumeric);
    p.addParameter('yTicks', [], @isnumeric);
    p.addParameter('titleLabel', '', @ischar);
    p.parse(varargin{:});
    
    [~,idx] = max(fittedActivationMap(:));
    [midRow, midCol] = ind2sub(size(fittedActivationMap),idx);
    activationSlice = squeeze(fittedActivationMap(midRow,:));
    offset = visualizeSpaceRange(1);
    gain =  (visualizeSpaceRange(2)-visualizeSpaceRange(1));
    
    cmap = parula(1024);
    hold(ax, 'on');
    faceAlpha = 1;
    edgeColor = [0.5 0.5 0.5];
    imagesc(ax, thePSFsupportDegs, thePSFsupportDegs, zeros(numel(thePSFsupportDegs), numel(thePSFsupportDegs)));
    plotConeApertures(ax, conePosDegs, coneAperturesDegs, integratedConeBasedSubgegionActivation, edgeColor, faceAlpha);
    [X,Y] = meshgrid(hiResPSFsupportDegs, hiResPSFsupportDegs);
    contour(ax,X,Y, fittedActivationMap, 0.1:0.2:1.0, 'LineColor', [1 1 0], 'LineWidth', 1.5);
    line(ax, hiResPSFsupportDegs, offset+gain*activationSlice, 'Color', [1 1 0], 'LineWidth', 1.5);
    
    axis(ax, 'xy'); axis(ax, 'square'); box(ax, 'on');
    set(ax, 'XLim', visualizeSpaceRange, 'YLim', visualizeSpaceRange, 'CLim', [0 1], 'ZLim', [0 1]);
    
    set(ax, 'XTick', p.Results.xTicks, 'YTick', p.Results.yTicks);
    if (~isempty(p.Results.titleLabel))
        title(ax, sprintf('%s\n(integrated)', p.Results.titleLabel));
    end
    colormap(ax, cmap);
end


function plotConeApertures(ax, conePosDegs, coneAperturesDegs, coneActivations, edgeColor, faceAlpha)
    for iCone = 1:size(conePosDegs,1)
        xr = conePosDegs(iCone,1) + coneAperturesDegs(iCone)*0.5*cosd(0:20:360);
        yr = conePosDegs(iCone,2) + coneAperturesDegs(iCone)*0.5*sind(0:20:360);
        %patch(ax, xr, yr, cmap(floor(coneActivations(iCone)*1023)+1,:), 'EdgeColor', edgeColor, 'FaceAlpha', faceAlpha);
        col = coneActivations(iCone);
        f = 1:numel(xr);
        v = [xr(:) yr(:)];
        patch(ax, 'Faces', f, 'Vertices', v, 'FaceVertexCData', col, 'FaceColor', 'flat', 'FaceAlpha', faceAlpha, 'EdgeColor', edgeColor);
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

function [thePSF, thePSFsupportDegs] = generatePSFForPatchAtEccentricity(PolansSubjectID, ...
    imposedRefractionErrorDiopters, wavelengthSampling, micronsPerDegree, patchEccDegs)

    pupilDiamMM = 3.0;
    theOI = CronerKaplanRGCModel.generatePolansOpticsForDeconcolution(...
        PolansSubjectID, imposedRefractionErrorDiopters, pupilDiamMM, ...
        wavelengthSampling, micronsPerDegree, patchEccDegs, ...
        'eccentricityUnits', 'degrees');
    
    % Get PSF slice at target wavelength
    theOptics = oiGet(theOI, 'optics');
    targetWavelength = 550;
    thePSF = opticsGet(theOptics,'psf data',targetWavelength);

    % Extract support
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

