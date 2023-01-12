
function [theRMSEvector, theRotatedRF, theRFprofile, ...
          theVisualSTF, theSpatialFrequencySupport, ...
          theFittedSTFsurroundToCenterRcRatio, ...
          theFittedSTFsurroundToCenterIntegratedSensitivityRatio, ratioWeights] = performCronerKaplanSimulation(...
            theVisualRF, theRetinalRFcenterConeMap, theRetinalRFsurroundConeMap, ...
            currentRetinalPoolingParams, retinalConePoolingParams, ... 
            bestHorizontalResolutionRotationDegs, spatialSupportDegs, ...
            targetSTFmatchMode, theTargetSTF, targetRsRcRatio, targetIntSensSCRatio, ...
            visualRFcenterCharacteristicRadiusDegs, multiStartsNumDoGFit, figNo, figureName)

    % Rotate theVisualRFmap according to the rotation
    % that maximizes horizontal resolution of the targetVisualRFmap
    theRotatedRF = RetinaToVisualFieldTransformer.bestHorizontalResolutionRFmap(...
        theVisualRF, bestHorizontalResolutionRotationDegs);

    % Integrate along Y to generate the visual RF profile
    theRFprofile = sum(theRotatedRF,1);

    % Compute the visual STF corresponding to the visual RF profile
    [theSpatialFrequencySupport, theVisualSTF] = ...
        RetinaToVisualFieldTransformer.spatialTransferFunction(...
                spatialSupportDegs, theRFprofile);

    if (isempty(theTargetSTF))
        theRMSEvector = [];
        theFittedSTFsurroundToCenterRcRatio = [];
        theFittedSTFsurroundToCenterIntegratedSensitivityRatio = [];
    else

        % Compute RMSEvector
        switch  (targetSTFmatchMode) 
            case 'STFcurve'
               % RMSE based on the difference of the STF curves
               theRMSEvector = ((theVisualSTF(:) - theTargetSTF(:))).^2;
    
               theFittedSTFsurroundToCenterRcRatio = [];
               theFittedSTFsurroundToCenterIntegratedSensitivityRatio = [];
    
               figure(1555); clf;
               ax = subplot(1,3,1);
               plot(ax,theSpatialFrequencySupport, theVisualSTF, 'ko-', 'MarkerSize', 12, 'LineWidth', 1.0); hold on;
               plot(ax,theSpatialFrequencySupport, theTargetSTF, 'b-', 'LineWidth', 1.5);
               set(ax, 'XScale', 'log', 'XLim', [0.1 100]);
               grid(ax,'on')
               

            case 'STFDoGparams'
                % Fit DoG model to the visual STF
                [theFittedSTFDoGparams, theFittedVisualSTF] = RetinaToVisualFieldTransformer.fitDoGmodelToMeasuredSTF(...
                      theSpatialFrequencySupport, ...
                      theVisualSTF, ...
                      visualRFcenterCharacteristicRadiusDegs, ...
                      multiStartsNumDoGFit, ...
                      'rangeForRc', visualRFcenterCharacteristicRadiusDegs*[0.5 1 1.5]);

                theFittedSTFsurroundToCenterRcRatio = theFittedSTFDoGparams.finalValues(3);
                theFittedSTFsurroundToCenterIntegratedSensitivityRatio = theFittedSTFDoGparams.finalValues(2) * (theFittedSTFDoGparams.finalValues(3))^2;
   

                % Assemble the RMSEvector
                %theRMSEvector = ((theVisualSTF(:) - theTargetSTF(:))).^2;

                ratioWeights = [1 10 3];
                shapeResidual = mean(theVisualSTF(:) - theTargetSTF(:)) / mean(theTargetSTF(:));
        
                %correlationMatrix = corrcoef(theVisualSTF(:), theTargetSTF(:));
                %theRMSEvector = (1-correlationMatrix(1,2))^2;

                
                RsRcRatioResidual = theFittedSTFsurroundToCenterRcRatio/targetRsRcRatio - 1;
                SCintSensRatioResidual = theFittedSTFsurroundToCenterIntegratedSensitivityRatio/targetIntSensSCRatio - 1;
                
               
                theRMSEvector(1,1) = ratioWeights(1) * shapeResidual^2;
                theRMSEvector(2,1) = ratioWeights(2) * RsRcRatioResidual^2;
                theRMSEvector(3,1) = ratioWeights(3) * SCintSensRatioResidual^2;

                rfColorMap = brewermap(1024, 'greys');

                hFig = figure(figNo); clf;
                set(hFig, 'Position', [10 10 1580 760], 'Name', figureName);

                if (~isempty(theRetinalRFcenterConeMap))
                    ax = subplot(2,4,1);
                    xMid = round(size(theRetinalRFcenterConeMap,2)/2);
                    yMid = round(size(theRetinalRFcenterConeMap,1)/2);
                    xx = find(abs(spatialSupportDegs)< 0.1);
                    yy = xx;
                    theVisualizedRFportion = theRetinalRFcenterConeMap(yy,xx);
                    xCenterProfile = sum(theVisualizedRFportion,1);
                    yCenterProfile = sum(theVisualizedRFportion,2);
                    maxRFamplitude = max(theVisualizedRFportion(:));
                    imagesc(ax, spatialSupportDegs(xx), spatialSupportDegs(yy), theVisualizedRFportion);
                    set(gca, 'CLim', [0 maxRFamplitude], 'FontSize', 16);
                    axis(ax,'image')
                    set(ax, 'XTick', [-0.5:0.05:0.5], 'YTick', [-0.5:0.05:0.5]);
                    title(ax,'retinal RF center')
                    ylabel(ax, 'space (degs)');
                    xlabel(ax, 'space (degs)');
                end
                colormap(ax, rfColorMap);

                ax = subplot(2,4,2);
                if (~isempty(theRetinalRFsurroundConeMap))
                    theVisualizedRFportion = theRetinalRFsurroundConeMap(yy,xx);
                    xSurroundProfile = sum(theVisualizedRFportion,1);
                    ySurroundProfile = sum(theVisualizedRFportion,2);
                    imagesc(ax, spatialSupportDegs(xx), spatialSupportDegs(yy), theVisualizedRFportion);
                    set(gca, 'CLim', [0 0.05*maxRFamplitude], 'FontSize', 16);
                    axis(ax,'image')
                    set(ax, 'XTick', [-0.5:0.05:0.5], 'YTick', [-0.5:0.05:0.5]);
                    title(ax,'retinal RF surround');
                    xlabel(ax, 'space (degs)');
                end
                colormap(ax, rfColorMap);


                ax = subplot(2,4,3);
                scalingFactor = max(theTargetSTF)/max(theFittedVisualSTF.compositeSTF);
                
                pMeasured = plot(ax,theSpatialFrequencySupport, theVisualSTF, 'ko', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 12, 'LineWidth', 1.0); hold on;
                pMeasuredDoGfit = plot(ax,theSpatialFrequencySupport, theFittedVisualSTF.compositeSTF, 'k-', 'LineWidth', 3.0);
                
                pTarget = plot(ax,theSpatialFrequencySupport, theTargetSTF, 'r-', 'LineWidth', 1.5);
                plot(ax,theSpatialFrequencySupport, theFittedVisualSTF.compositeSTF*scalingFactor, 'k--', 'LineWidth', 3.0);
                
                legend(ax, [pMeasured pMeasuredDoGfit pTarget], {'current STF', 'DoG fit', 'target STF'}, 'Location', 'SouthWest');
                axis(ax, 'square')
                set(ax, 'XScale', 'log', 'XLim', [0.1 100], 'XTick', [0.1 0.3 1 3 10 30 100], ...
                    'XTickLabel', {'0.1', '0.3', '1', '3', '10', '30', '100'}, ...
                    'YLim', [0 1.1*max(theTargetSTF)], 'FontSize', 16);
                grid(ax, 'on')
                xtickangle(ax, 0);
                title(ax,sprintf('Rs/Rc: achieved=%2.2f, target=%2.2f\nintS/C: achieved=%2.2f, target=%2.2f)', ...
                    theFittedSTFsurroundToCenterRcRatio, targetRsRcRatio, ...
                    theFittedSTFsurroundToCenterIntegratedSensitivityRatio, ...
                    targetIntSensSCRatio), 'FontSize', 16);
                xlabel(ax, 'spatial frequency (c/deg)');


                % The 1D-profiles
                ax = subplot(2,4, 5);
                xp = plot(ax,spatialSupportDegs(xx), xCenterProfile, 'r-', 'LineWidth', 1.0); hold(ax, 'on');
                plot(ax,spatialSupportDegs(xx),-xSurroundProfile, 'r-', 'Color', [0.5 0 0], 'LineWidth', 1.0);
                yp = plot(ax,spatialSupportDegs(yy),yCenterProfile, 'b-','LineWidth', 1.0);
                plot(ax,spatialSupportDegs(yy),-ySurroundProfile, 'b-', 'Color', [0 0 0.5],'LineWidth', 1.0);
                set(ax, 'XTick', [-0.5:0.05:0.5])
                grid(ax, 'on');
                axis(ax, 'square');
                set(ax,'FontSize', 16);
                title(ax, 'x- (red) and y- (blue) profiles');
                xlabel(ax, 'space (degs)');

                ax = subplot(2,4, 6);
                RetinaToVisualFieldTransformer.visualizeRetinalSurroundModelParametersAndRanges(ax, theFittedSTFDoGparams);
                title('DoG model params')

                ax = subplot(2,4, [7 8]);
                theCurrentRetinalConePoolingParams = retinalConePoolingParams;
                theCurrentRetinalConePoolingParams.finalValues = currentRetinalPoolingParams;
                RetinaToVisualFieldTransformer.visualizeRetinalSurroundModelParametersAndRanges(ax, theCurrentRetinalConePoolingParams );
                title('surround cone pooling params')

                drawnow;
    
            otherwise
                    error('modelConstants.targetSTFmatchMode must be set to either ''STFcurve'' or ''STFDoGparams''.');
        end % switch
    end

end
