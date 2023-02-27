function visualizeConeMosaicSTFresponses(mosaicFileName, responsesFileName, varargin)
    % Parse input
    p = inputParser;
    p.addParameter('opticsParams', [], @(x)(isempty(x)||isstruct(x))); 
    p.addParameter('targetConePositions', [], @(x)(isempty(x)||(size(x,2) == 2)));
    p.parse(varargin{:});
   
    targetConePositions = p.Results.targetConePositions;
    opticsParams = p.Results.opticsParams;

    load(mosaicFileName, 'theMidgetRGCMosaic')

    if (~isempty(opticsParams))       % Generate the native optics
        theMidgetRGCMosaic.generateNativeOptics(opticsParams);
        visualizedWavelength = 550;
        micronsPerDegree = theMidgetRGCMosaic.inputConeMosaic.micronsPerDegree;
        thePSFdata = generateNativeOpticsPSFdata(...
            theMidgetRGCMosaic.theNativeOptics, ...
            visualizedWavelength, micronsPerDegree);
    else
        thePSFdata = [];
    end


    load(responsesFileName, 'theNativeOpticsParams', ...
        'theConeMosaicSTFresponses', 'theConeMosaicNullResponses', ...
         'orientationsTested', 'spatialFrequenciesTested', ...
         'spatialPhasesDegs', 'coneContrasts');

    coneIndicesWithZeroNullResponse = find(theConeMosaicNullResponses== 0);
    normalizingResponses = 1./theConeMosaicNullResponses;
    normalizingResponses(coneIndicesWithZeroNullResponse) = 0;
    normalizingResponses = reshape(normalizingResponses, [1 1 numel(normalizingResponses)]);

    hFig = figure(10);
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 1280 1280]);
    ax = subplot('Position', [0.03 0.04 0.97 0.94]);
    cMap = gray(1024);

    for iCone = 1:size(targetConePositions,1)
        theSingleConeAxes{iCone} = axes('Position', [0.12 + (iCone-1)*0.16 0.78 0.15 0.15]);
    end

    if (~isempty(targetConePositions))
        for iCone = 1:size(targetConePositions,1)
            d = bsxfun(@minus, theMidgetRGCMosaic.inputConeMosaic.coneRFpositionsDegs, targetConePositions(iCone,:));
            d = sum(d.^2,2);
            [~, theIdentifiedConeIndices(iCone)] = min(d(:));
        end
    else
        theIdentifiedConeIndices = [];
    end

    
    videoOBJ = VideoWriter('ConeMosaicSTFresponses.mp4', 'MPEG-4');
    videoOBJ.FrameRate = 10;
    videoOBJ.Quality = 100;
    videoOBJ.open();

    
    orientationIndicesToVisualize = [7];
    sfIndicesToVisualize = [3 5 7 9 11 13 14];

    % Go through all stimulus orientations
    for iOri = 7:7
        orientationDegs = orientationsTested(iOri);

        theConeMosaicContrastResponsesAtThisOrientation = ...
            bsxfun(@times, bsxfun(@minus, squeeze(theConeMosaicSTFresponses(iOri,:,:,:)), theConeMosaicNullResponses), ...
                 normalizingResponses);

        for iFreq = 1:numel(spatialFrequenciesTested)

            for iFrame = 1:numel(spatialPhasesDegs)

                for iCone = 1:size(targetConePositions,1)
        
                    switch (theMidgetRGCMosaic.inputConeMosaic.coneTypes(theIdentifiedConeIndices(iCone)))
                        case cMosaic.LCONE_ID
                            coneColor = theMidgetRGCMosaic.inputConeMosaic.lConeColor;
                        case cMosaic.MCONE_ID
                            coneColor = theMidgetRGCMosaic.inputConeMosaic.mConeColor;
                        case cMosaic.SCONE_ID
                            coneColor = theMidgetRGCMosaic.inputConeMosaic.sConeColor;
                    end


                    m = 1.05*max(max(max(abs(theConeMosaicContrastResponsesAtThisOrientation(iFreq,:,theIdentifiedConeIndices)))));
                    if ( m < 0.1)
                        yTickIncrement = 0.02;
                    elseif (m < 0.5)
                        yTickIncrement = 0.1;
                    else
                        yTickIncrement = 0.2;
                    end

                    theSingleConeResponse = squeeze(theConeMosaicContrastResponsesAtThisOrientation(iFreq,:, theIdentifiedConeIndices(iCone)));
                    plot(theSingleConeAxes{iCone}, 1:iFrame, theSingleConeResponse(1:iFrame), 'wo-', 'Color', coneColor, 'LineWidth', 1.0, 'MarkerSize', 8, 'MarkerFaceColor', coneColor);
                    set(theSingleConeAxes{iCone}, 'FontSize', 16, 'LineWidth', 1.0, 'XTick', [], 'YTick', -1:yTickIncrement:1, 'YTickLabel', sprintf('%2.2f\n', -1:yTickIncrement:1), 'XTickLabel', {}, 'TickDir', 'both', 'XLim',[0 numel(spatialPhasesDegs)+1],  'XColor', 'none', 'YColor' ,[1 1.0 0.5], 'Color', [0.1 0.1 0.1], 'YLim', m*[-1 1]);
                    if (iCone>1)
                        set(theSingleConeAxes{iCone},'YTickLabel', {});
                    else
                        ylabel(theSingleConeAxes{iCone}, 'modulation');
                    end
                    grid(theSingleConeAxes{iCone}, 'on')
                    box(theSingleConeAxes{iCone}, 'off');

                end

                theMidgetRGCMosaic.inputConeMosaic.visualize(...
                            'figureHandle', hFig, ...
                            'axesHandle',ax, ...
                            'activation', theConeMosaicContrastResponsesAtThisOrientation(iFreq,iFrame,:), ...
                            'activationColorMap', cMap, ...
                            'activationRange', [-1 1]*m, ...
                            'labelConesWithIndices', theIdentifiedConeIndices, ...
                            'backgroundColor', [0 0 0], ...
                            'domainVisualizationTicks', struct(...
                            'x', 0:0.2:5, 'y', -1.6:0.2:1.6), ...
                            'plotTitle', sprintf('%2.1f c/deg', spatialFrequenciesTested(iFreq)), ...
                            'plotTitleColor', [1 1 0.5]);
                 

                 thePSFdataMinus = thePSFdata;
                 thePSFdataMinus.supportXdegs = thePSFdataMinus.supportXdegs-1.6;
                 thePSFaxis1 = axes('Position', [0.42-0.25 0.08 0.2 0.20]);
                 theMidgetRGCMosaic.inputConeMosaic.visualize(...
                            'figureHandle', hFig, ...
                            'axesHandle',thePSFaxis1, ...
                            'visualizedConeAperture', 'lightCollectingAreaCharacteristicDiameter', ...
                            'withSuperimposedPSF', thePSFdataMinus, ...
                            'domainVisualizationLimits', [2.4-1.6 2.6-1.6 -0.1 0.1], ...
                            'domainVisualizationTicks', struct(...
                            'x', 0:0.2:5, 'y', -1.6:0.2:1.6), ...
                            'noXLabel', true, 'noYLabel', true, ...
                            'fontSize', 1, ...
                            'backgroundColor', [0 0 0]);

                 thePSFaxis2 = axes('Position', [0.42 0.08 0.2 0.20]);
                 theMidgetRGCMosaic.inputConeMosaic.visualize(...
                            'figureHandle', hFig, ...
                            'axesHandle',thePSFaxis2, ...
                            'visualizedConeAperture', 'lightCollectingAreaCharacteristicDiameter', ...
                            'withSuperimposedPSF', thePSFdata, ...
                            'domainVisualizationLimits', [2.4 2.6 -0.1 0.1], ...
                            'domainVisualizationTicks', struct(...
                            'x', 0:0.2:5, 'y', -1.6:0.2:1.6), ...
                            'noXLabel', true, 'noYLabel', true, ...
                            'fontSize', 1, ...
                            'backgroundColor', [0 0 0]);
                 

                 thePSFdataPlus = thePSFdata;
                 thePSFdataPlus.supportXdegs = thePSFdataPlus.supportXdegs+1.6;
                 thePSFaxis3 = axes('Position', [0.42+0.25 0.08 0.2 0.20]);
                 theMidgetRGCMosaic.inputConeMosaic.visualize(...
                            'figureHandle', hFig, ...
                            'axesHandle',thePSFaxis3, ...
                            'visualizedConeAperture', 'lightCollectingAreaCharacteristicDiameter', ...
                            'withSuperimposedPSF', thePSFdataPlus, ...
                            'domainVisualizationLimits', [2.4+1.6 2.6+1.6 -0.1 0.1], ...
                            'domainVisualizationTicks', struct(...
                            'x', 0:0.2:5, 'y', -1.6:0.2:1.6), ...
                            'noXLabel', true, 'noYLabel', true, ...
                            'fontSize', 1, ...
                            'backgroundColor', [0 0 0]);
                 
                 set(hFig, 'Color', [0 0 0]);
                 set(ax, 'XColor', [0.7 0.7 0.7], 'YColor', [0.7 0.7 0.7]);


                 drawnow;
                 videoOBJ.writeVideo(getframe(hFig));
            end

        end
    end
    videoOBJ.close();

end


function thePSFData = generateNativeOpticsPSFdata(theOI, visualizedWavelength, micronsPerDegree)

    optics = oiGet(theOI, 'optics');
    waves = opticsGet(optics, 'wave');

    psfSupportMicrons = opticsGet(optics,'psf support','um');
    xGridDegs = psfSupportMicrons{1}/micronsPerDegree;
    yGridDegs = psfSupportMicrons{2}/micronsPerDegree;

    thePSFData.supportXdegs = xGridDegs(1,:);
    thePSFData.supportYdegs = yGridDegs(:,1);

    
    [~,idx] = min(abs(visualizedWavelength-waves));
    targetWavelength = waves(idx);
    thePSFData.data = opticsGet(optics,'psf data',targetWavelength );
    
end

