function [modelConstants, retinalConePoolingParams, visualRcDegs] = computeOptimizationComponents(...
    obj, theRGCindex, visualizeComponents, exportedFittingProgressFolder)

    % Compute the visual and anatomical Rc for this RGC
    [visualRcDegs, anatomicalRcDegs, indicesOfConesPooledByTheRFcenter, weightsOfConesPooledByTheRFcenter] = ...
        computeRFcenterVisualRc(obj, theRGCindex, visualizeComponents, exportedFittingProgressFolder);

    % Retinal cone pooling model constants
    modelConstants = struct();
    retinalConePoolingParams = struct();

    modelConstants.rmseWeightForRsRcResidual = obj.rmseWeightForRsRcResidual;
    modelConstants.rmseWeightForSCintSensResidual = obj.rmseWeightForSCintSensResidual;

    % The cone mosaic and the spectrally-weighted PSFs
    %modelConstants.theConeMosaic = obj.theRGCMosaic.inputConeMosaic;

    % The connectable cone types to the center and surroud
    modelConstants.surroundConnectableConeTypes = obj.retinalRFmodelParams.surroundConnectableConeTypes;
    modelConstants.centerConnectableConeTypes = obj.retinalRFmodelParams.centerConnectableConeTypes;
    modelConstants.inputConeMosaicConeTypes = obj.theRGCMosaic.inputConeMosaic.coneTypes;

    % Maximum support for the surround, in degrees, takes as twice the C&K
    % surrounds at the cell's eccentricity
    radialEccentricityForThisRGC = sqrt(sum(obj.theRGCMosaic.rgcRFpositionsDegs(theRGCindex,:).^2,2));
    
    % Surround cones up to 2 * surround characteristic radius of
    % Croner&Kaplan at this eccentricity
    modelConstants.maxSurroundSupportDegs = MosaicPoolingOptimizer.maxSurroundSupportFactor  * ...
                RGCmodels.CronerKaplan.constants.surroundCharacteristicRadiusFromFitToPandMcells(radialEccentricityForThisRGC);

    % Compute the RF center position
    RFcenterPos = mean(obj.theRGCMosaic.inputConeMosaic.coneRFpositionsDegs(indicesOfConesPooledByTheRFcenter,:),1);

    % Compute the distances of ALL cones in the input cone mosaic from the
    % RF center - THIS IS CACHED HERE
    coneDistancesFromRFCenterSquared = sum(bsxfun(@minus, obj.theRGCMosaic.inputConeMosaic.coneRFpositionsDegs, RFcenterPos).^2,2);
    coneDistancesFromRFCenter = sqrt(coneDistancesFromRFCenterSquared);
    modelConstants.cachedData.surroundConeIndices = find(coneDistancesFromRFCenter <= modelConstants.maxSurroundSupportDegs);
    modelConstants.cachedData.coneDistancesFromRFCenter = coneDistancesFromRFCenter(modelConstants.cachedData.surroundConeIndices);


    switch (obj.retinalRFmodelParams.conePoolingModel)
        case {
                MosaicPoolingOptimizer.ArbitraryCenter_DoubleExpH1cellIndex1Surround, ...
                MosaicPoolingOptimizer.ArbitraryCenter_DoubleExpH1cellIndex2Surround, ...
                MosaicPoolingOptimizer.ArbitraryCenter_DoubleExpH1cellIndex3Surround, ...
                MosaicPoolingOptimizer.ArbitraryCenter_DoubleExpH1cellIndex4Surround, ...
                }

            % From the 4 cells in Figure 6 of Packer & Dacey (2002)
            RnarrowToRwideRatios  = MosaicPoolingOptimizer.PackerDacey2002_H1params.RnarrowToRwideRatios;
            NWvolumeRatios        = MosaicPoolingOptimizer.PackerDacey2002_H1params.NWvolumeRatios;

            modelConstants.retinalConePoolingModel = obj.retinalRFmodelParams.conePoolingModel;
            modelConstants.weightsComputeFunctionHandle = @MosaicPoolingOptimizer.conePoolingCoefficientsForArbitraryCenterDoubleExpSurround;
            modelConstants.indicesOfCenterCones = indicesOfConesPooledByTheRFcenter;
            modelConstants.weightsOfCenterCones = weightsOfConesPooledByTheRFcenter;

            
            % Range for RwDegs, based to characteristic radius of cones and # of cones in RF center
            RwDegsInitial    = 4.00 * 1/mean(RnarrowToRwideRatios) * anatomicalRcDegs * sqrt(2.3) * sqrt(numel(modelConstants.indicesOfCenterCones));
            RwDegsLowerBound = max([0.02 0.01*RwDegsInitial]);
            RwDegsUpperBound = min([modelConstants.maxSurroundSupportDegs 4*RwDegsInitial]);

            %                                        Kc      Ks/KcRatio    narrowToWideFieldVolumeRatio  RwideDegs            RnarrowToRwideRatio
            retinalConePoolingParams.names =         {'Kc',  'KsKcRatio',  'VnVwRatio',                  'RwDegs',             'RnRwRatio'};
            retinalConePoolingParams.scaling =       {'log', 'log',        'log',                        'linear',                'log'};
            retinalConePoolingParams.initialValues = [1.       0.06        mean(NWvolumeRatios)           RwDegsInitial         mean(RnarrowToRwideRatios)];
            retinalConePoolingParams.lowerBounds   = [0.5      0.005       min(NWvolumeRatios)            RwDegsLowerBound      min(RnarrowToRwideRatios)];
            retinalConePoolingParams.upperBounds   = [2        1e0         max(NWvolumeRatios)            RwDegsUpperBound      max(RnarrowToRwideRatios)];

            if (contains(modelConstants.retinalConePoolingModel, 'H1cellIndex1'))
                H1cellIndex = 1;
                fprintf('Fitting using params from the FIRST H1 cell of Dacey and Packer\n');
            elseif (contains(modelConstants.retinalConePoolingModel, 'H1cellIndex2'))
                H1cellIndex = 2;
                fprintf('Fitting using params from the SECOND H1 cell of Dacey and Packer\n');
            elseif (contains(modelConstants.retinalConePoolingModel, 'H1cellIndex3'))
                H1cellIndex = 3;
                fprintf('Fitting using params from the THIRD H1 cell of Dacey and Packer\n');
            elseif (contains(modelConstants.retinalConePoolingModel, 'H1cellIndex4'))
                H1cellIndex = 4;
                fprintf('Fitting using params from the FOURTH H1 cell of Dacey and Packer\n');
            else
                error('Could not determine H1cellIndex from the retinal cone pooling model: ''%s''.\n', ...
                    modelConstants.retinalConePoolingModel);
            end

            idx = find(ismember(retinalConePoolingParams.names, 'VnVwRatio'));
            measuredValue = NWvolumeRatios(H1cellIndex);
            retinalConePoolingParams.initialValues(idx) = measuredValue;
            retinalConePoolingParams.lowerBounds(idx) = measuredValue * (1 - obj.retinalRFmodelParams.H1parameterTolerance);
            retinalConePoolingParams.upperBounds(idx) = measuredValue * (1 + obj.retinalRFmodelParams.H1parameterTolerance);

            idx = find(ismember(retinalConePoolingParams.names, 'RnRwRatio'));
            measuredValue = RnarrowToRwideRatios(H1cellIndex);
            retinalConePoolingParams.initialValues(idx) = measuredValue;
            retinalConePoolingParams.lowerBounds(idx) = measuredValue * (1 - obj.retinalRFmodelParams.H1parameterTolerance);
            retinalConePoolingParams.upperBounds(idx) = measuredValue * (1 + obj.retinalRFmodelParams.H1parameterTolerance);
           

        otherwise
            error('Retinal cone pooling model ''%s'' is not implemented\n', obj.retinalRFmodelParams.conePoolingModel)
    end % switch (obj.retinalRFmodelParams.retinalConePoolingModel)
end

function [visualRcDegs, anatomicalRcDegs, inputConeIndices, inputConeWeights] = ...
    computeRFcenterVisualRc(obj, theTargetRGCindex, visualizeComponents, exportedFittingProgressFolder)

    inputConeIndices = find(squeeze(obj.theRGCMosaic.rgcRFcenterConeConnectivityMatrix(:,theTargetRGCindex))>0.001);
    inputConeWeights = full(obj.theRGCMosaic.rgcRFcenterConeConnectivityMatrix(inputConeIndices, theTargetRGCindex));

    theCenterResponsesAcrossAllOrientationsAndSpatialFrequencies = sum(bsxfun(@times, ...
        obj.inputConeMosaicVisualSTFdata.responseModulations(:,:,:,inputConeIndices), ...
        reshape(inputConeWeights, [1 1 1 numel(inputConeIndices)])),4);

    % The STF of center
    orientationsTested = obj.inputConeMosaicVisualSTFdata.orientationsTested;
    spatialFrequenciesTested = obj.inputConeMosaicVisualSTFdata.spatialFrequenciesTested;

    [theOptimalCenterSTF, theCenterSTFsAcrossAllOrientations, theOptimalOrientation] = ...
        MosaicPoolingOptimizer.optimalSTFfromResponsesToAllOrientationsAndSpatialFrequencies( ...
                    orientationsTested, spatialFrequenciesTested, ...
                    theCenterResponsesAcrossAllOrientationsAndSpatialFrequencies);

    if (visualizeComponents)

        domainVisualizationCenter = obj.theRGCMosaic.rgcRFpositionsDegs(theTargetRGCindex,:);
        domainVisualizationLimits = [domainVisualizationCenter(1) domainVisualizationCenter(1) domainVisualizationCenter(2) domainVisualizationCenter(2)] + ...
                                    0.12*[-1 1 -1 1];
        domainVisualizationTicks = struct('x', -10:0.05:10, 'y', -10:0.05:10);

        hFig = figure(996); clf;
        ff = MSreadyPlot.figureFormat('1x1 small');
        theAxes = MSreadyPlot.generateAxes(hFig,ff);
        set(hFig, 'Color', [1 1 1]);

        % Generate and set the optics
        obj.theRGCMosaic.setTheOptics(obj.inputConeMosaicVisualSTFdata.metaData.theNativeOpticsParams);

        visualizedWavelength = 550;
        micronsPerDegree = obj.theRGCMosaic.inputConeMosaic.micronsPerDegree;
        thePSFdata = mRGCMosaic.generateOpticsPSFdataForVisualization(...
            obj.theRGCMosaic.theNativeOptics, ...
            visualizedWavelength, micronsPerDegree);

        % Center the visualized PSF at the [ domainVisualizationLimits(1), domainVisualizationLimits(3)]
        thePSFdata.supportXdegs = thePSFdata.supportXdegs - obj.theRGCMosaic.inputConeMosaic.eccentricityDegs(1) + ...
            domainVisualizationCenter(1);
        thePSFdata.supportYdegs = thePSFdata.supportYdegs - obj.theRGCMosaic.inputConeMosaic.eccentricityDegs(2) + ...
            domainVisualizationCenter(2);

        % Show just the PSF
        obj.theRGCMosaic.inputConeMosaic.visualize(...
            'figureHandle', hFig,...
            'axesHandle', theAxes{1,1}, ...
            'domain', 'degrees', ...
            'domainVisualizationLimits', domainVisualizationLimits, ...
            'domainVisualizationTicks', domainVisualizationTicks, ...
            'backgroundColor', [1 1 1], ...
            'visualizeCones', ~true,...
            'visualizedConeAperture', 'lightCollectingAreaCharacteristicDiameter', ...
            'withSuperimposedPSF', thePSFdata, ...
            'fontAngle', ff.axisFontAngle, ...
            'fontSize', ff.fontSize);
        
        if (~isempty(exportedFittingProgressFolder))
            rawFiguresRoot = exportedFittingProgressFolder;
            pdfFileName = fullfile(rawFiguresRoot,'RGCmosaicInputConeMosaicAtSurroundOptimizationNode.pdf');
            NicePlot.exportFigToPDF(pdfFileName, hFig, 300);
        end

        hFig = figure(997); clf;
        ff = MSreadyPlot.figureFormat('1x1 small');
        theAxes = MSreadyPlot.generateAxes(hFig,ff);
        set(hFig, 'Color', [1 1 1]);

        obj.theRGCMosaic.visualize(...
            'figureHandle', hFig,...
            'axesHandle', theAxes{1,1}, ...
            'domainVisualizationLimits', domainVisualizationLimits, ...
            'domainVisualizationTicks', domainVisualizationTicks, ...
            'identifyInputCones', true, ...
            'identifiedInputConeIndices', inputConeIndices, ...
            'identifyPooledCones', ~true, ...
            'identifiedConeAperture', 'lightCollectingArea4sigma', ...
            'plottedRFoutlineLineWidth', 1.5, ...
            'centerSubregionContourSamples', 30, ...
            'labelRGCsWithIndices', theTargetRGCindex, ...
            'labeledRGCsLineWidth', ff.lineWidth, ...
            'labeledRGCsColor', [1 0 0], ...
            'backgroundColor', [1 1 1], ...
            'fontAngle', ff.axisFontAngle, ...
            'fontSize', ff.fontSize);

        if (~isempty(exportedFittingProgressFolder))
            rawFiguresRoot = exportedFittingProgressFolder;
            pdfFileName = fullfile(rawFiguresRoot,'RGCmosaicAtSurroundOptimizationNode.pdf');
            NicePlot.exportFigToPDF(pdfFileName, hFig, 300);
        end




        % Generate 2D visual STF
        SF = [];
        ORI = [];
        STF = [];

        for iOri = 1:numel(orientationsTested)
            for iSF = 1:numel(spatialFrequenciesTested)
                SF = cat(1, SF, spatialFrequenciesTested(iSF));
                ORI = cat(1, ORI, orientationsTested(iOri));
                STF = cat(1, STF, theCenterSTFsAcrossAllOrientations(iOri, iSF));
                if (orientationsTested(iOri) == 0)
                    SF = cat(1, SF, spatialFrequenciesTested(iSF));
                    ORI = cat(1, ORI, 180);
                    STF = cat(1, STF, theCenterSTFsAcrossAllOrientations(iOri, iSF));
                end
            end
        end

        sfX = SF .* cosd(ORI);
        sfY = SF .* sind(ORI);

        F = scatteredInterpolant(sfX, sfY, STF, 'natural');

        hiResSFx = -70:0.5:70;
        hiResSFy = (0:0.5:70)'; 
        [hiResSFxx, hiResSFyy] = meshgrid(hiResSFx, hiResSFy);
        theCenterSTFsAcrossAllOrientationsHiRes = F(hiResSFxx(:),hiResSFyy(:));
        theCenterSTFsAcrossAllOrientationsHiRes = theCenterSTFsAcrossAllOrientationsHiRes / max(theCenterSTFsAcrossAllOrientationsHiRes(:));
        theCenterSTFsAcrossAllOrientationsHiRes = reshape(theCenterSTFsAcrossAllOrientationsHiRes, [numel(hiResSFy) numel(hiResSFx)]);
        theCenterSTFsAcrossAllOrientationsHiRes(theCenterSTFsAcrossAllOrientationsHiRes<0.001) = 0.001;

        % Make it symmetric
        theCenterSTFsAcrossAllOrientationsHiRes = cat(1, ...
            flipud(fliplr(theCenterSTFsAcrossAllOrientationsHiRes)), ...
            theCenterSTFsAcrossAllOrientationsHiRes(2:end,:) ...
            );

        hFig = figure(999); clf;
        ff = MSreadyPlot.figureFormat('1x1 small');
        theAxes = MSreadyPlot.generateAxes(hFig,ff);
        set(hFig, 'Color', [1 1 1]);

        contourLevelsLogMag = log10([1 2 4 8 15 30 55 95]/100);
        cLut = brewermap(100, 'OrRd');
        cLut(1,:) = cLut(1,:)*0.6 + 0.4*[1 1 1];

        [M,c] = contourf(theAxes{1,1}, hiResSFx, hiResSFx, log10(theCenterSTFsAcrossAllOrientationsHiRes), ...
            contourLevelsLogMag);
        c.LineWidth = ff.axisLineWidth*1.5;
        c.LineStyle = '-';
        c.EdgeColor = [0.3 0.3 0.5];

        % 20 % contour in red
        hold(theAxes{1,1}, 'on');
        contourLevels2 = log10(MosaicPoolingOptimizer.highSFAttenuationFactorForOptimalOrientation)*[1 1];
        [M2,c2] = contourf(theAxes{1,1}, hiResSFx, hiResSFx, log10(theCenterSTFsAcrossAllOrientationsHiRes), ...
            contourLevels2);
        c2.LineStyle = '-';
        c2.LineWidth = ff.axisLineWidth*2;
        c2.EdgeColor = [0 0 0];
        c2.FaceAlpha = 0;

        axis(theAxes{1,1}, 'image');
        axis(theAxes{1,1}, 'xy');
        set(theAxes{1,1}, 'TickDir', 'both', 'XTick', -100:20:100, 'YTick', -100:20:100);
        set(theAxes{1,1}, 'XLim', [-70 70], 'YLim', [-70 70], 'CLim', [min(contourLevelsLogMag) max(contourLevelsLogMag)]);

        % Font size
        set(theAxes{1,1}, 'FontSize', ff.fontSize, 'FontAngle', ff.axisFontAngle);

        % axis color and width
        set(theAxes{1,1}, 'XColor', ff.axisColor, 'YColor', ff.axisColor, 'Color', [1 1 1], 'LineWidth', ff.axisLineWidth);
        grid(theAxes{1,1}, 'on')
        box(theAxes{1,1}, 'off')

        xtickangle(theAxes{1,1}, 0);

        xlabel(theAxes{1,1}, 'spatial frequency, x (c/deg)')
        ylabel(theAxes{1,1}, 'spatial frequency, y (c/deg)');

        
        colormap(theAxes{1,1}, cLut);
        
        colorbar('Ticks',contourLevelsLogMag,...
                 'TickLabels',sprintf('%.0f%%\n', 100*(10.^contourLevelsLogMag)), ...
                 'FontSize', ff.fontSize-2);

        if (~isempty(exportedFittingProgressFolder))
            rawFiguresRoot = exportedFittingProgressFolder;
            pdfFileName = fullfile(rawFiguresRoot,'RFcenter2DVisualSTF.pdf');
            NicePlot.exportFigToPDF(pdfFileName, hFig, 300);

            % Generate paper-ready figures (scaled versions of the figures in
            % the rawFiguresRoot directory) which are stored in the PaperReady folder
            dd = strrep(rawFiguresRoot, 'Raw', 'cpdf');

            PLOSdirectory = '/Users/nicolas/Documents/4_LaTeX/PLOS2023-Overleaf/matlabFigureCode';
            commandString = sprintf('%s -args %s/generatePLOSOnePaperReadyFigures.txt', dd, PLOSdirectory);
            system(commandString);

            pause
        end

    end


    % An estimate of the anatomical RcDegs
    anatomicalRcDegs = sqrt(numel(inputConeWeights)) * ...
                       obj.theRGCMosaic.inputConeMosaic.coneApertureToConeCharacteristicRadiusConversionFactor * ...
                       mean(obj.theRGCMosaic.inputConeMosaic.coneApertureDiametersDegs(inputConeIndices));

    % Some initial estimate of the visual Rc (this is clearly off as it
    % does not account the effect of optics), but it is better than
    % nothing. The Gaussian STF model is only a 2-parameter model (gain,
    % Rc), so it should be able to find the Rc without getting stuct to a
    % local minimum
    initialRcDegs = anatomicalRcDegs;
    
    % Fit the visual STF with a DoG model
    fittedParamsStruct = MosaicPoolingOptimizer.fitGaussianToSubregionSTF(...
                      obj.inputConeMosaicVisualSTFdata.spatialFrequenciesTested, ...
                      theOptimalCenterSTF, ...
                      initialRcDegs, ...
                      [], ...
                      obj.multiStartsNumDoGFit);

    visualRcDegs = fittedParamsStruct.finalValues(2);
end