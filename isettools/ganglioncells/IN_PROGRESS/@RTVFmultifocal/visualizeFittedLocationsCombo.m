function visualizeFittedLocationsCombo(mosaicDirectory, figNo, theMidgetRGCmosaic, ...
    theRTFVTobjList, theOpticsPositionGrid, theConesNumPooledByTheRFcenterGrid, ...
    varargin)

    validComponents = {...
        'retinal quantal efficiencies', ...
        'PSFs', ...
        'retinal L-center RF subregions', ...
        'retinal M-center RF subregions', ...
        'retinal L-center RF subregions & L-weighted PSF', ...
        'retinal M-center RF subregions & M-weighted PSF', ...
        'L-center STFs & L-weighted PSF', ...
        'M-center STFs & M-weighted PSF', ...
        'L-center composite RFs', ...
        'M-center composite RFs' ...
    };

    p = inputParser;
    p.addParameter('topPanelInfo', 'retinal quantal efficiencies', @(x)(ismember(x, validComponents)));
    p.addParameter('bottomPanelInfo', 'PSFs', @(x)(ismember(x, validComponents)));
    p.addParameter('visualizedSpatialRangeArcMin', 10, @isscalar);
    p.addParameter('figPostfix', '', @ischar);

    p.parse(varargin{:});
    topPanelInfo = p.Results.topPanelInfo;
    bottomPanelInfo = p.Results.bottomPanelInfo;
    figPostfix = p.Results.figPostfix;
    visualizedSpatialRangeArcMin = p.Results.visualizedSpatialRangeArcMin;

    visualizationRequiresComputationOfRFs = false;
    if (contains(topPanelInfo, 'RF')) || (contains(bottomPanelInfo, 'RF')) || ...
       (contains(topPanelInfo, 'STF')) || (contains(bottomPanelInfo, 'STF')) 
        visualizationRequiresComputationOfRFs = true;
    end

    % Find how many different #of center cones are fitted
    conesNumPooled = unique(theConesNumPooledByTheRFcenterGrid);
    xCoords = squeeze(theMidgetRGCmosaic.inputConeMosaic.coneRFpositionsDegs(:,1));
    yCoords = squeeze(theMidgetRGCmosaic.inputConeMosaic.coneRFpositionsDegs(:,2));
    xLimsDegs = [min(xCoords) max(xCoords)];
    yLimsDegs = [min(yCoords) max(yCoords)];

    psfRangeDegs = visualizedSpatialRangeArcMin/60;
    rfRangeDegs = psfRangeDegs;
 
    for iConesNumPooled = 1:numel(conesNumPooled)

        theConesNumPooled = conesNumPooled(iConesNumPooled);
        [eccentricitySamplingGrid, theRTVFobjIndicesForThisGrid] = RTVFmultifocal.subGridSpatialCoordsForConesNumPooled(theConesNumPooled, ...
            theConesNumPooledByTheRFcenterGrid, theOpticsPositionGrid);

        % Determine which RGCs have this many center cones
        rgcIndices = find(sum(theMidgetRGCmosaic.rgcRFcenterConePoolingMatrix,1) == theConesNumPooled);

        for iSpatialPosition = 1:size(eccentricitySamplingGrid,1)
            theOpticsPosition = eccentricitySamplingGrid(iSpatialPosition,:);

            % Extract the corresponding RTVF object index
            dd = bsxfun(@minus,eccentricitySamplingGrid, theOpticsPosition);
            [~, idx] = min(sum(dd.^2,2));
            theRTVFobjIndex = theRTVFobjIndicesForThisGrid(idx);
            theRTVFobj = theRTFVTobjList{theRTVFobjIndex};

            if (visualizationRequiresComputationOfRFs)
                % Compute the fitted L- and M-cone RF data
                fittedLconeModelRFdata = extractRFandSTFdata(theRTVFobj, theRTVFobj.LconeRFcomputeStruct);
                fittedMconeModelRFdata = extractRFandSTFdata(theRTVFobj, theRTVFobj.MconeRFcomputeStruct);
            end

            % Compute the retinal L- and M-cone quantal efficiencies at this
            % optical position, taking into account the input cone mosaic's
            % variation in MP density with eccentricity
            [LquantalEfficiency, MquantalEffiency,wavelengthSupport] = ...
                RTVF.LMconeSpectralWeightings(theRTVFobj.coneMosaic, theOpticsPosition);

            % Extract the PSFdata at this spatial position
            thePSFData.psfSupportXdegs = theRTVFobj.spectrallyWeightedPSFData.psfSupportXdegs;
            thePSFData.psfSupportYdegs = theRTVFobj.spectrallyWeightedPSFData.psfSupportYdegs;
            thePSFData.LconeWeighted = ...
                conv2(theRTVFobj.spectrallyWeightedPSFData.LconeWeighted, ...
                      theRTVFobj.coneApertureBlurKernel, 'same');
            thePSFData.MconeWeighted = ...
                conv2(theRTVFobj.spectrallyWeightedPSFData.MconeWeighted, ...
                      theRTVFobj.coneApertureBlurKernel, 'same');
            
            hFig = figure(figNo+iConesNumPooled*100 + iSpatialPosition); clf;
            figName = sprintf('RTVF_%d_%s.pdf', theRTVFobjIndex, figPostfix);
            ff = MSreadyPlot.figureFormat('1x2 large');
            set(hFig, 'Position', [10 10 ff.figureSize(1) ff.figureSize(2)], ...
                'Color', [1 1 1], 'Name', figName);
            
            % Left panel: RGC positions and RTVF sampling points (red)
            ax = subplot('Position', ff.subplotPosVectors(1,1).v);
            MSreadyPlot.renderRFpositions(ax, theMidgetRGCmosaic.rgcRFpositionsDegs(rgcIndices,:), ...
                xLimsDegs, yLimsDegs, sprintf('%d-cone RF centers', theConesNumPooled), ff);
    
            hold(ax, 'on');
            % All positions in red
            plot(ax, eccentricitySamplingGrid(:,1), eccentricitySamplingGrid(:,2), '+', ...
               'Color', [1 0.5 0.5], 'LineWidth', ff.lineWidth*2, 'MarkerSize', ff.markerSize);
            plot(ax, eccentricitySamplingGrid(:,1), eccentricitySamplingGrid(:,2), '+', ...
               'Color', [1 0 0], 'LineWidth', ff.lineWidth, 'MarkerSize', ff.markerSize);
     
            % The current position in blue
            plot(ax, theOpticsPosition(1), theOpticsPosition(2), '+', ...
               'Color', [0.5 1 1], 'LineWidth', ff.lineWidth*4, 'MarkerSize', 2*ff.markerSize);
            plot(ax, theOpticsPosition(1), theOpticsPosition(2), '+', ...
               'Color', [0 0 1], 'LineWidth', ff.lineWidth*2 , 'MarkerSize', 2*ff.markerSize);
    
            % Top panel plots
            switch (topPanelInfo)
                case 'retinal quantal efficiencies'
                    ax = subplot('Position',  [0.62 0.6 0.3 0.35]);
                    MSreadyPlot.renderConeFundamentals(ax, wavelengthSupport, ...
                        LquantalEfficiency, MquantalEffiency, [], 'retinal quantal efficiency (qe)', ff);
                
                case 'L-center STFs & L-weighted PSF'
                    % Top left panel: computed and DoG model fit of the visual STF (L) at current location
                    ax = subplot('Position',  [0.52 0.62 0.14 0.36]);

                    MSreadyPlot.renderSTF(ax, ...
                        fittedLconeModelRFdata.theSTFdata.spatialFrequencySupport, ...
                        fittedLconeModelRFdata.theSTFdata.visualSTF, ...
                        fittedLconeModelRFdata.theSTFdata.fittedDoGModelToVisualSTF.compositeSTF, ...
                        fittedLconeModelRFdata.theSTFdata.fittedDoGModelToVisualSTF.centerSTF, ...
                        fittedLconeModelRFdata.theSTFdata.fittedDoGModelToVisualSTF.surroundSTF, ...
                        '', ...
                        {'achieved STF', 'fitted DoG STF', 'fitted center STF', 'fitted surround STF'}, ff, ...
                        'noYLabel', true);

                    % Top middle panel: correspondence between achieved and desired DoG ratios at current location
                    targetSCintSensRatio = theRTVFobj.targetVisualRFDoGparams.surroundToCenterIntegratedSensitivityRatio;
                    targetRsRcRatio = theRTVFobj.targetVisualRFDoGparams.surroundToCenterRcRatio;
                    achievedRsRcRatio = fittedLconeModelRFdata.theSTFdata.fittedDoGModelRsRcRatio;
                    achievedSCintSensRatio = fittedLconeModelRFdata.theSTFdata.fittedDoGModelSCIntSensRatio;
                    
                    ax = subplot('Position', [0.71 0.62 0.10 0.29]);
                    MSreadyPlot.renderPerformance(ax, ...
                        targetRsRcRatio, targetSCintSensRatio, ...
                        achievedRsRcRatio, achievedSCintSensRatio, ...
                        ff);

                    % Top right panel: L-cone qe weighted PSF at current location
                    ax = subplot('Position', [0.84 0.62 0.14 0.29]);
                    MSreadyPlot.render2DPSF(ax, thePSFData.psfSupportXdegs, thePSFData.psfSupportYdegs, ...
                        thePSFData.LconeWeighted, psfRangeDegs, sprintf('L-cone q.e.-\nweighted PSF'), ff, ...
                        'noXLabel', false, 'noYLabel', true);

                case 'retinal L-center RF subregions & L-weighted PSF'
                    % Top left panel: RF center (L) at current location
                    ax = subplot('Position',  [0.52 0.62 0.14 0.29]);
                    MSreadyPlot.renderRFsubregion(ax, ...
                        fittedLconeModelRFdata.xSupportDegs, fittedLconeModelRFdata.ySupportDegs, ...
                        fittedLconeModelRFdata.theRetinalRFcenterConeMap, ...
                        rfRangeDegs, sprintf('retinal center\n(L-center)'), ff, ...
                        'withLineWeightingFunction', fittedLconeModelRFdata.theRetinalRFcenterLineWeightingFunction);

                    
                    % Top middle panel: RF surround (L) at current location
                    ax = subplot('Position', [0.68 0.62 0.14 0.29]);
                    MSreadyPlot.renderRFsubregion(ax, ...
                        fittedLconeModelRFdata.xSupportDegs, fittedLconeModelRFdata.ySupportDegs, ...
                        fittedLconeModelRFdata.theRetinalRFsurroundConeMap, ...
                        rfRangeDegs, sprintf('retinal surround\n(L-center)'), ff, ...
                        'withLineWeightingFunction', fittedLconeModelRFdata.theRetinalRFsurroundLineWeightingFunction, ...
                        'noYLabel', true);

                    % Top right panel: L-cone qe weighted PSF at current location
                    ax = subplot('Position', [0.84 0.62 0.14 0.29]);
                    MSreadyPlot.render2DPSF(ax, thePSFData.psfSupportXdegs, thePSFData.psfSupportYdegs, ...
                        thePSFData.LconeWeighted, psfRangeDegs, sprintf('L-cone q.e.-\nweighted PSF'), ff, ...
                        'noXLabel', false, 'noYLabel', true);

                case 'retinal L-center RF subregions'
                    % Top left panel: RF center (L) at current location
                    ax = subplot('Position',  [0.55 0.59 0.22 0.37]);
                    MSreadyPlot.renderRFsubregion(ax, ...
                        fittedLconeModelRFdata.xSupportDegs, fittedLconeModelRFdata.ySupportDegs, ...
                        fittedLconeModelRFdata.theRetinalRFcenterConeMap, ...
                        rfRangeDegs, 'retinal RF center (L-center)', ff, ...
                        'withLineWeightingFunction', fittedLconeModelRFdata.theRetinalRFcenterLineWeightingFunction);

                    % Top right panel: RF surround (L) at current location
                    ax = subplot('Position', [0.78 0.59 0.22 0.37]);
                    MSreadyPlot.renderRFsubregion(ax, ...
                        fittedLconeModelRFdata.xSupportDegs, fittedLconeModelRFdata.ySupportDegs, ...
                        fittedLconeModelRFdata.theRetinalRFsurroundConeMap, ...
                        rfRangeDegs, 'retinal RF surround (L-center)', ff, ...
                        'withLineWeightingFunction', fittedLconeModelRFdata.theRetinalRFsurroundLineWeightingFunction, ...
                        'noYLabel', true);


                case 'retinal M-center RF subregions'
                    % Top left panel: retinal RF center (M) at current location
                    ax = subplot('Position',  [0.55 0.59 0.22 0.37]);
                    MSreadyPlot.renderRFsubregion(ax, ...
                        fittedMconeModelRFdata.xSupportDegs, fittedMconeModelRFdata.ySupportDegs, ...
                        fittedMconeModelRFdata.theRetinalRFcenterConeMap, ...
                        rfRangeDegs, 'retinal RF center (M-center)', ff, ...
                        'withLineWeightingFunction', fittedMconeModelRFdata.theRetinalRFcenterLineWeightingFunction);

                    % Top right panel: retinal RF surround (M) at current location
                    ax = subplot('Position', [0.78 0.59 0.22 0.37]);
                    MSreadyPlot.renderRFsubregion(ax, ...
                        fittedMconeModelRFdata.xSupportDegs, fittedMconeModelRFdata.ySupportDegs, ...
                        fittedMconeModelRFdata.theRetinalRFsurroundConeMap, ...
                        rfRangeDegs, 'retinal RF surround (M-center)', ff, ...
                        'withLineWeightingFunction', fittedMconeModelRFdata.theRetinalRFsurroundLineWeightingFunction, ...
                        'noYLabel', true);


                case 'PSFs'
                    % Top left panel: L-cone PSF at current location
                    ax = subplot('Position',  [0.55 0.59 0.22 0.37]);
                    MSreadyPlot.render2DPSF(ax, thePSFData.psfSupportXdegs, thePSFData.psfSupportYdegs, ...
                        thePSFData.LconeWeighted, psfRangeDegs, 'L-cone qe - weighted PSF', ff, ...
                        'noXLabel', false, 'noYLabel', false);
                    
                    % Top right panel: M-cone PSF at current location
                    ax = subplot('Position', [0.78 0.59 0.22 0.37]);
                    MSreadyPlot.render2DPSF(ax, thePSFData.psfSupportXdegs, thePSFData.psfSupportYdegs, ...
                        thePSFData.MconeWeighted, psfRangeDegs, 'M-cone qe - weighted PSF', ff, ...
                        'noXLabel', false, 'noYLabel', true);


                case 'L-center composite RFs'
                    % Top left panel: retinal composite RF (M) at current location
                    ax = subplot('Position',  [0.55 0.59 0.22 0.37]);
                    MSreadyPlot.renderRFcompositeMap(ax, ...
                        fittedLconeModelRFdata.xSupportDegs, fittedLconeModelRFdata.ySupportDegs, ...
                        fittedLconeModelRFdata.theRetinalRFcenterConeMap - ...
                        fittedLconeModelRFdata.theRetinalRFsurroundConeMap, ...
                        rfRangeDegs, 'retinal composite RF (L-center)', ff);

                    % Top right panel: visual composite RF (M) at current location
                    ax = subplot('Position', [0.78 0.59 0.22 0.37]);
                    MSreadyPlot.renderRFcompositeMap(ax, ...
                        fittedLconeModelRFdata.xSupportDegs, fittedLconeModelRFdata.ySupportDegs, ...
                        fittedLconeModelRFdata.theVisualRFmap, ...
                        rfRangeDegs, 'visual RF (L-center)', ff);

                case 'M-center composite RFs'
                    % Top left panel: retinal composite RF (M) at current location
                    ax = subplot('Position',  [0.55 0.59 0.22 0.37]);
                    MSreadyPlot.renderRFcompositeMap(ax, ...
                        fittedMconeModelRFdata.xSupportDegs, fittedMconeModelRFdata.ySupportDegs, ...
                        fittedMconeModelRFdata.theRetinalRFcenterConeMap - ...
                        fittedMconeModelRFdata.theRetinalRFsurroundConeMap, ...
                        rfRangeDegs, 'retinal composite RF (M-center)', ff);

                    % Top right panel: visual composite RF (M) at current location
                    ax = subplot('Position', [0.78 0.07 0.22 0.37]);
                    MSreadyPlot.renderRFcompositeMap(ax, ...
                        fittedMconeModelRFdata.xSupportDegs, fittedMconeModelRFdata.ySupportDegs, ...
                        fittedMconeModelRFdata.theVisualRFmap, ...
                        rfRangeDegs, 'visual RF (M-center)', ff);
                    
                otherwise
                    error('Top panel info ''%s'' not implented', topPanelInfo);
            end

            % Bottom panel plots
            switch (bottomPanelInfo)
                case 'retinal quantal efficiencies'
                    ax = subplot('Position',  [0.62 0.08 0.3 0.35]);
                    MSreadyPlot.renderConeFundamentals(ax, wavelengthSupport, ...
                        LquantalEfficiency, MquantalEffiency, [], 'retinal quantal efficiency (qe)', ff);
                
                case 'M-center STFs & M-weighted PSF'
                    % Top left panel: computed and DoG model fit of the visual STF (M) at current location
                    ax = subplot('Position',  [0.52 0.14 0.14 0.36]);

                    MSreadyPlot.renderSTF(ax, ...
                        fittedMconeModelRFdata.theSTFdata.spatialFrequencySupport, ...
                        fittedMconeModelRFdata.theSTFdata.visualSTF, ...
                        fittedMconeModelRFdata.theSTFdata.fittedDoGModelToVisualSTF.compositeSTF, ...
                        fittedMconeModelRFdata.theSTFdata.fittedDoGModelToVisualSTF.centerSTF, ...
                        fittedMconeModelRFdata.theSTFdata.fittedDoGModelToVisualSTF.surroundSTF, ...
                        '', ...
                        {'achieved STF', 'fitted DoG STF', 'fitted center STF', 'fitted surround STF'}, ff, ...
                        'noYLabel', true);

                    % Top middle panel: correspondence between achieved and desired DoG ratios at current location
                    targetSCintSensRatio = theRTVFobj.targetVisualRFDoGparams.surroundToCenterIntegratedSensitivityRatio;
                    targetRsRcRatio = theRTVFobj.targetVisualRFDoGparams.surroundToCenterRcRatio;
                    
                    achievedRsRcRatio = fittedMconeModelRFdata.theSTFdata.fittedDoGModelRsRcRatio;
                    achievedSCintSensRatio = fittedMconeModelRFdata.theSTFdata.fittedDoGModelSCIntSensRatio;
                    
                    ax = subplot('Position', [0.71 0.14 0.10 0.29]);
                    MSreadyPlot.renderPerformance(ax, ...
                        targetRsRcRatio, targetSCintSensRatio, ...
                        achievedRsRcRatio, achievedSCintSensRatio, ...
                        ff);

                    % Top right panel: M-cone qe weighted PSF at current location
                    ax = subplot('Position', [0.84 0.14 0.14 0.29]);
                    MSreadyPlot.render2DPSF(ax, thePSFData.psfSupportXdegs, thePSFData.psfSupportYdegs, ...
                        thePSFData.MconeWeighted, psfRangeDegs, sprintf('M-cone q.e.-\nweighted PSF'), ff, ...
                        'noXLabel', false, 'noYLabel', true);


                case 'retinal M-center RF subregions & M-weighted PSF'
                    % Bottom left panel: RF center (M) at current location
                    ax = subplot('Position',  [0.52 0.14 0.14 0.29]);
                    MSreadyPlot.renderRFsubregion(ax, ...
                        fittedMconeModelRFdata.xSupportDegs, fittedMconeModelRFdata.ySupportDegs, ...
                        fittedMconeModelRFdata.theRetinalRFcenterConeMap, ...
                        rfRangeDegs, sprintf('retinal center\n(M-center)'), ff, ...
                        'withLineWeightingFunction', fittedMconeModelRFdata.theRetinalRFcenterLineWeightingFunction);

                    
                    % Bottom middle panel: RF surround (M) at current location
                    ax = subplot('Position', [0.68 0.14 0.14 0.29]);
                    MSreadyPlot.renderRFsubregion(ax, ...
                        fittedMconeModelRFdata.xSupportDegs, fittedMconeModelRFdata.ySupportDegs, ...
                        fittedMconeModelRFdata.theRetinalRFsurroundConeMap, ...
                        rfRangeDegs, sprintf('retinal surround\n(M-center)'), ff, ...
                        'withLineWeightingFunction', fittedMconeModelRFdata.theRetinalRFsurroundLineWeightingFunction, ...
                        'noYLabel', true);

                    % Bottom right panel: M-cone qe weighted PSF at current location
                    ax = subplot('Position', [0.84 0.14 0.14 0.29]);
                    MSreadyPlot.render2DPSF(ax, thePSFData.psfSupportXdegs, thePSFData.psfSupportYdegs, ...
                        thePSFData.MconeWeighted, psfRangeDegs, sprintf('M-cone q.e.-\nweighted PSF'), ff, ...
                        'noXLabel', false, 'noYLabel', true);

                case 'retinal L-center RF  subregions'
                    % Bottom left panel: RF center (L) at current location
                    ax = subplot('Position',  [0.55 0.07 0.22 0.37]);
                    MSreadyPlot.renderRFsubregion(ax, ...
                        fittedLconeModelRFdata.xSupportDegs, fittedLconeModelRFdata.ySupportDegs, ...
                        fittedLconeModelRFdata.theRetinalRFcenterConeMap, ...
                        rfRangeDegs, 'retinal RF center (L-center)', ff, ...
                        'withLineWeightingFunction', fittedLconeModelRFdata.theRetinalRFcenterLineWeightingFunction);

                    % Bottom right panel: RF surround (L) at current location
                    ax = subplot('Position', [0.78 0.07 0.22 0.37]);
                    MSreadyPlot.renderRFsubregion(ax, ...
                        fittedLconeModelRFdata.xSupportDegs, fittedLconeModelRFdata.ySupportDegs, ...
                        fittedLconeModelRFdata.theRetinalRFsurroundConeMap, ...
                        rfRangeDegs, 'retinal RF surround (L-center)', ff, ...
                        'withLineWeightingFunction', fittedLconeModelRFdata.theRetinalRFsurroundLineWeightingFunction, ...
                        'noYLabel', true);


                case 'retinal M-center RF subregions'
                    % Bottom left panel: retinal RF center (M) at current location
                    ax = subplot('Position',  [0.55 0.07 0.22 0.37]);
                    MSreadyPlot.renderRFsubregion(ax, ...
                        fittedMconeModelRFdata.xSupportDegs, fittedMconeModelRFdata.ySupportDegs, ...
                        fittedMconeModelRFdata.theRetinalRFcenterConeMap, ...
                        rfRangeDegs, 'retinal RF center (M-center)', ff, ...
                        'withLineWeightingFunction', fittedMconeModelRFdata.theRetinalRFcenterLineWeightingFunction);

                    % Bottom right panel: retinal RF surround (M) at current location
                    ax = subplot('Position', [0.78 0.07 0.22 0.37]);
                    MSreadyPlot.renderRFsubregion(ax, ...
                        fittedMconeModelRFdata.xSupportDegs, fittedMconeModelRFdata.ySupportDegs, ...
                        fittedMconeModelRFdata.theRetinalRFsurroundConeMap, ...
                        rfRangeDegs, 'retinal RF surround (M-center)', ff, ...
                        'withLineWeightingFunction', fittedMconeModelRFdata.theRetinalRFsurroundLineWeightingFunction, ...
                        'noYLabel', true);


                case 'PSFs'
                    % Bottom left panel: L-cone PSF at current location
                    ax = subplot('Position',  [0.55 0.07 0.22 0.37]);
                    MSreadyPlot.render2DPSF(ax, thePSFData.psfSupportXdegs, thePSFData.psfSupportYdegs, ...
                        thePSFData.LconeWeighted, psfRangeDegs, 'L-cone qe - weighted PSF', ff, ...
                        'noXLabel', false, 'noYLabel', false);

                    % Bottom right panel: M-cone PSF at current location
                    ax = subplot('Position', [0.78 0.07 0.22 0.37]);
                    MSreadyPlot.render2DPSF(ax, thePSFData.psfSupportXdegs, thePSFData.psfSupportYdegs, ...
                        thePSFData.MconeWeighted, psfRangeDegs, 'M-cone qe - weighted PSF', ff, ...
                        'noXLabel', false, 'noYLabel', true);

                case 'L-center composite RFs'
                    % Bottom left panel: retinal composite RF (M) at current location
                    ax = subplot('Position',  [0.55 0.07 0.22 0.37]);
                    MSreadyPlot.renderRFcompositeMap(ax, ...
                        fittedLconeModelRFdata.xSupportDegs, fittedLconeModelRFdata.ySupportDegs, ...
                        fittedLconeModelRFdata.theRetinalRFcenterConeMap - ...
                        fittedLconeModelRFdata.theRetinalRFsurroundConeMap, ...
                        rfRangeDegs, 'retinal composite RF (L-center)', ff);

                    % Bottom right panel: visual composite RF (M) at current location
                    ax = subplot('Position', [0.78 0.07 0.22 0.37]);
                    MSreadyPlot.renderRFcompositeMap(ax, ...
                        fittedLconeModelRFdata.xSupportDegs, fittedLconeModelRFdata.ySupportDegs, ...
                        fittedLconeModelRFdata.theVisualRFmap, ...
                        rfRangeDegs, 'visual RF (L-center)', ff);

                case 'M-center composite RFs'
                    % Bottom left panel: retinal composite RF (M) at current location
                    ax = subplot('Position',  [0.55 0.07 0.22 0.37]);
                    MSreadyPlot.renderRFcompositeMap(ax, ...
                        fittedMconeModelRFdata.xSupportDegs, fittedMconeModelRFdata.ySupportDegs, ...
                        fittedMconeModelRFdata.theRetinalRFcenterConeMap - ...
                        fittedMconeModelRFdata.theRetinalRFsurroundConeMap, ...
                        rfRangeDegs, 'retinal composite RF (M-center)', ff);

                    % Bottom right panel: visual composite RF (M) at current location
                    ax = subplot('Position', [0.78 0.07 0.22 0.37]);
                    MSreadyPlot.renderRFcompositeMap(ax, ...
                        fittedMconeModelRFdata.xSupportDegs, fittedMconeModelRFdata.ySupportDegs, ...
                        fittedMconeModelRFdata.theVisualRFmap, ...
                        rfRangeDegs, 'visual RF (M-center)', ff);

                otherwise
                    error('Bottom panel info ''%s'' not implented', bottomPanelInfo);
            end
            
            drawnow;
            NicePlot.exportFigToPDF(fullfile(mosaicDirectory,figName), hFig, 300);
        end % iSpatialPosition
    end % iConesNumPooled

end


function dataOut = extractRFandSTFdata(theRTVFTobj, theConeSpecificRFcomputeStruct)

    % Compute the visual RF map
    [theVisualRFmap, theRetinalRFcenterConeMap, ...
     theRetinalRFsurroundConeMap] = theRTVFTobj.visualRFfromRetinalConePooling(...
                theConeSpecificRFcomputeStruct.modelConstants, ...
                theConeSpecificRFcomputeStruct.retinalConePoolingParams.finalValues);
  
    % Compute the STF data for the visual RF map
    theSTFdata = theRTVFTobj.visualRFmapPropertiesFromCronerKaplanAnalysis(theVisualRFmap);

    dataOut = struct();
    dataOut.xSupportDegs = theConeSpecificRFcomputeStruct.modelConstants.spatialSupportDegs(:,1);
    dataOut.ySupportDegs = theConeSpecificRFcomputeStruct.modelConstants.spatialSupportDegs(:,2);
    dataOut.theRetinalRFcenterConeMap = theRetinalRFcenterConeMap;
    dataOut.theRetinalRFsurroundConeMap = theRetinalRFsurroundConeMap;
    dataOut.theVisualRFmap = theVisualRFmap;
    dataOut.theSTFdata = theSTFdata;

    centerLineWeightingFunction = sum(theRetinalRFcenterConeMap,1);
    surroundLineWeightingFunction = sum(theRetinalRFsurroundConeMap,1);
    maxLineWeightingFunction = max(centerLineWeightingFunction);
    dataOut.theRetinalRFcenterLineWeightingFunction = centerLineWeightingFunction / maxLineWeightingFunction;
    dataOut.theRetinalRFsurroundLineWeightingFunction = -surroundLineWeightingFunction / maxLineWeightingFunction;

end
