
function [theRMSEvector, theRotatedRF, theRFprofile, ...
          theVisualSTF, theSpatialFrequencySupport, ...
          theFittedSTFsurroundToCenterRcRatio, ...
          theFittedSTFsurroundToCenterIntegratedSensitivityRatio, ratioWeights] = performCronerKaplanSimulation(...
            theVisualRF, theRetinalRFcenterConeMap, theRetinalRFsurroundConeMap, ...
            currentRetinalPoolingParams, retinalConePoolingParams, ... 
            bestHorizontalResolutionRotationDegs, sfSupport, ...
            targetSTFmatchMode, theTargetSTF, targetRsRcRatio, targetIntSensSCRatio, ...
            visualRFcenterCharacteristicRadiusDegs, multiStartsNumDoGFit, figNo)

    % Rotate theVisualRFmap according to the rotation
    % that maximizes horizontal resolution of the targetVisualRFmap
    theRotatedRF = RetinaToVisualFieldTransformer.bestHorizontalResolutionRFmap(...
        theVisualRF, bestHorizontalResolutionRotationDegs);

    % Integrate along Y to generate the visual RF profile
    theRFprofile = sum(theRotatedRF,1);

    % Compute the visual STF corresponding to the visual RF profile
    [theSpatialFrequencySupport, theVisualSTF] = ...
        RetinaToVisualFieldTransformer.spatialTransferFunction(...
                sfSupport, theRFprofile);

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

                
                hFig = figure(figNo); clf;
                set(hFig, 'Position', [10 10 1000 1200]);
                if (~isempty(theRetinalRFcenterConeMap))
                    ax = subplot(3,2,1);
                    xMid = round(size(theRetinalRFcenterConeMap,2)/2);
                    yMid = round(size(theRetinalRFcenterConeMap,1)/2);
                    xx = xMid + (-(round(xMid/10)):1:(round(xMid/10)));
                    yy = yMid + (-(round(yMid/10)):1:(round(yMid/10)));
                    imagesc(ax, theRetinalRFcenterConeMap(yy,xx));
                    axis(ax,'image')
                    title(ax,'retinal RF center')
                end

                ax = subplot(3,2,2);
                if (~isempty(theRetinalRFsurroundConeMap))
                    imagesc(ax, theRetinalRFsurroundConeMap(yy,xx));
                    axis(ax,'image')
                    title(ax,'current retinal RF surround')
                end
                colormap(gray);


                ax = subplot(3,2,3);
                pMeasured = plot(ax,theSpatialFrequencySupport, theVisualSTF, 'ko', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 12, 'LineWidth', 1.0); hold on;
                pMeasuredDoGfit = plot(ax,theSpatialFrequencySupport, theFittedVisualSTF.compositeSTF, 'k-', 'LineWidth', 3.0);
                pTarget = plot(ax,theSpatialFrequencySupport, theTargetSTF, 'r-', 'LineWidth', 1.5);
                legend(ax, [pMeasured pMeasuredDoGfit pTarget], {'current', 'current (DoG fit)', 'target'}, 'Location', 'SouthWest');

                set(ax, 'XScale', 'log', 'XLim', [0.1 100]);
                grid(ax, 'on')
                title(ax,sprintf('Rs/Rc= (current/target: %2.2f/%2.2f)\nintS/C= (current/target: %2.2f/%2.2f)', ...
                    theFittedSTFsurroundToCenterRcRatio, targetRsRcRatio, ...
                    theFittedSTFsurroundToCenterIntegratedSensitivityRatio, ...
                    targetIntSensSCRatio), 'FontSize', 16);
                
                ax = subplot(3,2,5);
                RetinaToVisualFieldTransformer.visualizeRetinalSurroundModelParametersAndRanges(ax, theFittedSTFDoGparams);
                title('current DoG fit params')

                ax = subplot(3,2,6);
                theCurrentRetinalConePoolingParams = retinalConePoolingParams;
                theCurrentRetinalConePoolingParams .finalValues = currentRetinalPoolingParams;
                RetinaToVisualFieldTransformer.visualizeRetinalSurroundModelParametersAndRanges(ax, theCurrentRetinalConePoolingParams );
                title('current retinal cone pooling params')

                drawnow;
    
            otherwise
                    error('modelConstants.targetSTFmatchMode must be set to either ''STFcurve'' or ''STFDoGparams''.');
        end % switch
    end

end
