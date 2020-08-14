function deconvolutionStruct = performDeconvolutionAnalysis(conesNumInRFcenterTested, ...
    conePosDegs, coneAperturesDegs, thePSF, thePSFsupportDegs)
    
    % Deconvolution data is a container
    deconvolutionStruct.data = containers.Map();  
    
    % Generate the cone aperture profile
    coneApertureProfile = generateConeApertureProfile(thePSFsupportDegs, coneAperturesDegs);
    
    % For some numbers of cone inputs, e.g. 2-5, we examine how the PSF interacts with several possible combinations of 
    % these inputs. Below we generate these possible input combinations
    rfCenterPoolingSchemes = generatePoolingSchemes(conesNumInRFcenterTested, coneAperturesDegs, conePosDegs);
    
    inputConeIndices = 1:size(conePosDegs,1);
    coneMask = generateRetinalConeImage(inputConeIndices, conePosDegs, coneApertureProfile, thePSFsupportDegs);
    coneMask(coneMask>0.0001) = 1;
            
    % Whether to visualize the fits (for debuging purposes mainly)
    visualizeFits = ~true;
    
    for poolingSchemeIndex = 1:numel(rfCenterPoolingSchemes)
        
        inputConeIndicesForAllCombinations = rfCenterPoolingSchemes{poolingSchemeIndex}.inputConeIndices;
           
        for inputConeCombinationIndex = 1:size(inputConeIndicesForAllCombinations,1)
            fprintf('Analyzing %d of %d spatial input schemes for a %d-cone center subregion.\n',  ...
                inputConeCombinationIndex, size(inputConeIndicesForAllCombinations,1), rfCenterPoolingSchemes{poolingSchemeIndex}.coneInputsNum);
      
            % Retrieve the indices of cones feeding into the RF center
            inputConeIndices = inputConeIndicesForAllCombinations(inputConeCombinationIndex,:);
            
            % Compute the retinal cone image (all input cones lit-up)
            retinalConeImage = generateRetinalConeImage(inputConeIndices, conePosDegs, coneApertureProfile, thePSFsupportDegs);
            
            % Compute the visual cone image = retinal cone image * PSF
            visualConeImage = conv2(retinalConeImage, thePSF, 'same');
            
            % Integrate within cone apertures
            [retinalConeActivations, visualConeActivations, X, Y] = integrateWithinConeApertures(...
                retinalConeImage, visualConeImage, coneMask, conePosDegs, coneAperturesDegs, thePSFsupportDegs);

            % Normalize activations with respect to the retinal domain
            visualConeActivations = visualConeActivations / max(retinalConeActivations(:));
            retinalConeActivations = retinalConeActivations / max(retinalConeActivations(:));
            
            
            % Fit integraded cone activation maps
            [rfSigmasRetinal(inputConeCombinationIndex,:), rfGainRetinal(inputConeCombinationIndex),  ...
                hiResRetinalConeActivationMap, hiResPSFsupportDegs] = ...
                fitActivationMap(retinalConeActivations, conePosDegs, coneAperturesDegs, thePSFsupportDegs);
            
            [rfSigmasVisual(inputConeCombinationIndex,:), rfGainVisual(inputConeCombinationIndex),  ...
                hiResVisualConeActivationMap] = ...
                fitActivationMap(visualConeActivations, conePosDegs, coneAperturesDegs, thePSFsupportDegs);
            
            % Visualize fits
            if (visualizeFits)
                visualizeDeconvolutionFits(inputConeCombinationIndex, thePSFsupportDegs, thePSF, ...
                    conePosDegs, retinalConeImage, visualConeImage, ...
                    retinalConeActivations, visualConeActivations, ...
                    hiResPSFsupportDegs, hiResRetinalConeActivationMap, hiResVisualConeActivationMap, ...
                    rfSigmasRetinal(inputConeCombinationIndex,:), rfGainRetinal(inputConeCombinationIndex), ...
                    rfSigmasVisual(inputConeCombinationIndex,:), rfGainVisual(inputConeCombinationIndex));
                pause
            end
            
        end %inputConeCombination
        
        % Log data
        keyName = sprintf('%d-coneInput',rfCenterPoolingSchemes{poolingSchemeIndex}.coneInputsNum);
        deconvolutionStruct.data(keyName) = struct(...
            'retinalGain', mean(rfGainRetinal), ...
            'visualGain',  mean(rfGainVisual), ...
            'minRetinalSigma', mean(min(rfSigmasRetinal,[],2)), ...
            'maxRetinalSigma', mean(max(rfSigmasRetinal,[],2)), ...
            'minVisualSigma', mean(min(rfSigmasVisual,[],2)), ...
            'maxVisualSigma', mean(max(rfSigmasVisual,[],2)) ...
        );
    end % poolingSchemeIndex      
end


function visualizeDeconvolutionFits(inputConeCombinationIndex, thePSFsupportDegs, thePSF, ...
                conePosDegs, retinalConeImage, visualConeImage, ...
                retinalConeActivations, visualConeActivations, ...
                hiResPSFsupportDegs, hiResRetinalConeActivationMap, hiResVisualConeActivationMap, ...
                rfSigmasRetinal, rfGainRetinal, rfSigmasVisual, rfGainVisual)

    hFig = figure(inputConeCombinationIndex); clf;
    set(hFig, 'Position', [100 100 14 13], ...
        'Name', sprintf('Deconvolution fits for input layout #%d', inputConeCombinationIndex));

    theAxesGrid = plotlab.axesGrid(hFig, ...
            'leftMargin', 0.04, ...
            'bottomMargin', 0.04, ...
            'rightMargin', 0.01, ...
            'widthMargin', 0.05, ...
            'heightMargin', 0.08, ...
            'topMargin', 0.03, ...
            'rowsNum', 2, ...
            'colsNum', 2);

        
    ax = theAxesGrid{1,1};
    renderConeImage(ax,thePSFsupportDegs, retinalConeImage, thePSF, 'retinal cone image');

    ax = theAxesGrid{1,2};
    renderConeImage(ax,thePSFsupportDegs, visualConeImage, thePSF, 'visual cone image');

    ax = theAxesGrid{2,1};
    renderConeActivationImage(ax,retinalConeActivations, conePosDegs, ...
        hiResPSFsupportDegs, hiResRetinalConeActivationMap, ...
        sprintf('gain: %2.2f, sigmas: (%2.2f'',%2.2f'')', ...
        rfGainRetinal, rfSigmasRetinal(1)*60, rfSigmasRetinal(2)*60));

    ax = theAxesGrid{2,2};
    renderConeActivationImage(ax,visualConeActivations, conePosDegs, ...
        hiResPSFsupportDegs, hiResVisualConeActivationMap, ...
        sprintf('gain: %2.2f, sigmas: (%2.2f'',%2.2f'')', ...
        rfGainVisual, rfSigmasVisual(1)*60, rfSigmasVisual(2)*60));
    
end

function renderConeImage(ax,thePSFsupportDegs, coneImage, thePSF, figTitle)
    imagesc(ax,thePSFsupportDegs, thePSFsupportDegs, coneImage);
    hold(ax, 'on');
    contour(ax,thePSFsupportDegs, thePSFsupportDegs, thePSF/max(thePSF(:)), 0.2:0.2:1);
    plot(ax,[-1 1], [0 0], 'r-');
    plot(ax,[0 0],[-1 1], 'r-');
    set(ax, 'XLim', 0.1*[-1 1], 'YLim', 0.1*[-1 1], 'XTick', (-0.1:0.05:0.1), 'YTick', (-0.1:0.05:0.1));
    axis(ax, 'square'); axis(ax, 'xy');
    title(ax,figTitle);
end

function renderConeActivationImage(ax,retinalConeActivations, conePosDegs, ...
        hiResPSFsupportDegs, hiResRetinalConeActivationMap, figTitle)
    
    surf(ax,hiResPSFsupportDegs, hiResPSFsupportDegs, hiResRetinalConeActivationMap);
    hold(ax, 'on');
    stem3(ax,conePosDegs(:,1), conePosDegs(:,2),retinalConeActivations, 'filled');
    set(ax, 'XLim', 0.1*[-1 1], 'YLim', 0.1*[-1 1], 'XTick', (-0.1:0.05:0.1), 'YTick', (-0.1:0.05:0.1), 'CLim', [0 1], 'ZLim', [0 1]);
    axis(ax,'square');  axis(ax, 'xy');
    title(ax, figTitle);
end

function [rfSigmas, rfGain, hiResConeActivationMap, hiResPSFsupportDegs] = ...
    fitActivationMap(coneActivations, conePosDegs, coneAperturesDegs, thePSFsupportDegs)

    deltaX = thePSFsupportDegs(2)-thePSFsupportDegs(1);
    minConeCharacteristicRadiusDegs = 0.5*min(coneAperturesDegs)/3;
    
    % Fit Gaussian to cone activations
    functionName = 'elliptical Gaussian';
    [fittedParams, rfFunction] = ...
        fitElliptical2DGausianToRF(conePosDegs(:,1), conePosDegs(:,2), coneActivations(:), ...
        deltaX/40, minConeCharacteristicRadiusDegs, [0 0], functionName);
    
    
    % Retrieve fit params
    rfGain = fittedParams(1);
    if (strcmp(functionName, 'elliptical Gaussian'))
        rfSigmas = [fittedParams(4) fittedParams(5)];
    else
        rfSigmas = fittedParams(4)*[1 1];
    end
    
    % Generate a high-res fitted map
    maxPSFsupport = deltaX*100;
    hiResPSFsupportDegs = -maxPSFsupport:deltaX*2:maxPSFsupport;
    [XXX,YYY] = meshgrid(hiResPSFsupportDegs, hiResPSFsupportDegs);
    xyData(:,1) = XXX(:);
    xyData(:,2) = YYY(:);
    fittedActivationMap = rfFunction(fittedParams,xyData);
    hiResConeActivationMap = reshape(fittedActivationMap, [numel(hiResPSFsupportDegs) numel(hiResPSFsupportDegs)]);
end


function [retinalConeActivations, visualConeActivations, X, Y] = integrateWithinConeApertures(...
    retinalConeImage, visualConeImage, coneMask, conePosDegs, coneAperturesDegs, thePSFsupportDegs)

    conesNum = size(conePosDegs,1);
    retinalConeActivations = zeros(1, conesNum);
    visualConeActivations = zeros(1, conesNum);
    [X,Y] = meshgrid(thePSFsupportDegs,thePSFsupportDegs);
    
    retinalConeImage = retinalConeImage .* coneMask;
    visualConeImage = visualConeImage .* coneMask;
    
    for iCone = 1:conesNum
        xr = conePosDegs(iCone,1) + coneAperturesDegs(iCone)*0.5*cosd(0:2:360);
        yr = conePosDegs(iCone,2) + coneAperturesDegs(iCone)*0.5*sind(0:2:360);
        [in,on] = inpolygon(X(:), Y(:), xr, yr);
        withinAperturePixelIndices = find((in==true));
        if (numel(withinAperturePixelIndices) == 0)
            retinalConeActivations(iCone) = 0;
            visualConeActivations(iCone) = 0;
        else
            netRetinalConeActivation = sum(retinalConeImage(withinAperturePixelIndices));
            netVisualConeActivation  = sum(visualConeImage(withinAperturePixelIndices));
            retinalConeActivations(iCone) = netRetinalConeActivation;
            visualConeActivations(iCone) = netVisualConeActivation;
        end
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

function coneApertureProfile = generateConeApertureProfile(thePSFsupportDegs, coneAperturesDegs)
    
    % Spatial support for the cone aperture
    spatialSampleSize = thePSFsupportDegs(2)-thePSFsupportDegs(1);
    theConeApertureSupportDegs = (spatialSampleSize:spatialSampleSize:0.6*max(coneAperturesDegs));
    theConeApertureSupportDegs = [-fliplr(theConeApertureSupportDegs) 0 theConeApertureSupportDegs];
    [Xcone, Ycone] = meshgrid(theConeApertureSupportDegs, theConeApertureSupportDegs);
    
    % Cone aperture profile
    apertureShape = 'flattopGaussian';
    coneCharacteristicRadiusDegs = 0.5*mean(coneAperturesDegs)/3;
    if (strcmp(apertureShape, 'disk'))
        r = sqrt(Xcone.^2+Ycone.^2);
        coneApertureProfile = zeros(size(Xcone));
        coneApertureProfile(r<=0.5*mean(coneAperturesDegs)) = 1;
    elseif (strcmp(apertureShape, 'flattopGaussian'))
        coneApertureProfile = (exp(-(Xcone/coneCharacteristicRadiusDegs).^2) .* exp(-(Ycone/coneCharacteristicRadiusDegs).^2));
        flatTopThreshold = exp(-5)*exp(-5);
        coneApertureProfile(coneApertureProfile<flatTopThreshold) = 0;
        coneApertureProfile = coneApertureProfile .^ 0.01;
    else  
        coneApertureProfile = exp(-(Xcone/coneCharacteristicRadiusDegs).^2) .* exp(-(Ycone/coneCharacteristicRadiusDegs).^2);
    end
end

function poolingSchemes = generatePoolingSchemes(conesNumInRFcenterTested, coneAperturesDegs, conePosDegs)

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
    end
end