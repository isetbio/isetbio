function rmseSequence = displayFittingProgress(hFigProgress, videoOBJ, rmseSequence, ...
                RsRcRatioResidual, SCintSensRatioResidual, theCurrentRMSE, ...
                retinalConePoolingParams, currentRetinalPoolingParamValues, ...
                theCurrentSTFdata, ...
                spatialSupportDegsX, ...
                spatialSupportDegsY, ...
                theCurrentRetinalRFcenterConeMap, ...
                theCurrentRetinalRFsurroundConeMap)

            % Update the sequence of RMSEs
            rmseSequence(:,size(rmseSequence,2)+1) = abs([RsRcRatioResidual SCintSensRatioResidual theCurrentRMSE]);

            currentRetinalConePoolingParamsStruct = retinalConePoolingParams;
            currentRetinalConePoolingParamsStruct.finalValues = currentRetinalPoolingParamValues;

            figure(hFigProgress);
            ax = subplot(2,4,1);
            RTVF.visualizeFittedModelParametersAndRanges(ax, currentRetinalConePoolingParamsStruct);
    
            ax = subplot(2,4,5);
            RTVF.visualizeFittedModelParametersAndRanges(ax, theCurrentSTFdata.fittedDoGModelParams);

            ax = subplot(2,4,[2 6]);
            cla(ax);
            iterations = 1:size(rmseSequence,2);
            
            plot(ax,iterations, rmseSequence(1,:), 'r-', 'MarkerSize', 12, 'LineWidth', 1.0); hold(ax, 'on')
            plot(ax,iterations, rmseSequence(2,:), 'b-', 'MarkerSize', 12, 'LineWidth', 1.0); 
            plot(ax,iterations, rmseSequence(3,:), 'k-', 'MarkerSize', 12, 'LineWidth', 1.0); 
            legend(ax, {'Rs/Rc', 'int S/C', 'total'}, 'Location', 'NorthOutside', 'Orientation','horizontal');
            set(ax, 'YLim', max(abs(rmseSequence(:)))*[1e-5 1.1], 'YScale', 'log', 'FontSize', 14);
            xlabel(ax, 'iteration');
            ylabel(ax, 'RMSE');

            if (theCurrentRMSE == min(squeeze(rmseSequence(3,:)))) 
                if (~isempty(theCurrentRetinalRFcenterConeMap))
                    ax = subplot(2,4,3);
                    imagesc(ax, spatialSupportDegsX, spatialSupportDegsY, theCurrentRetinalRFcenterConeMap);
                    axis(ax, 'image');
                    title(sprintf('RF at iteration %d', size(rmseSequence,2)));
                    ax = subplot(2,4,7);
                    imagesc(ax, spatialSupportDegsX, spatialSupportDegsY, theCurrentRetinalRFsurroundConeMap);
                    axis(ax, 'image');
                    colormap(brewermap(1024, 'greys'));
    
                    ax = subplot(2,4,4);
                    cla(ax);
                    centerProfileX = sum(theCurrentRetinalRFcenterConeMap,1);
                    surroundProfileX = sum(theCurrentRetinalRFsurroundConeMap,1);
                    plotProfiles(ax, spatialSupportDegsX, centerProfileX, surroundProfileX);
    
                    centerProfileY = sum(theCurrentRetinalRFcenterConeMap,2);
                    surroundProfileY = sum(theCurrentRetinalRFsurroundConeMap,2);
                end


                ax = subplot(2,4,8);
                cla(ax);
                plot(ax,theCurrentSTFdata.spatialFrequencySupport, theCurrentSTFdata.visualSTF, 'ks', 'LineWidth', 1.5);
                hold(ax, 'on');
                plot(ax,theCurrentSTFdata.spatialFrequencySupport, theCurrentSTFdata.fittedDoGModelToVisualSTF.compositeSTF, 'k-', 'LineWidth', 1.5);
                plot(ax,theCurrentSTFdata.spatialFrequencySupport, theCurrentSTFdata.fittedDoGModelToVisualSTF.centerSTF, 'r-', 'LineWidth', 1.0);
                plot(ax,theCurrentSTFdata.spatialFrequencySupport, theCurrentSTFdata.fittedDoGModelToVisualSTF.surroundSTF, 'b-', 'LineWidth', 1.0);
                set(ax, 'XScale', 'log', 'XLim', [0.1 100], 'XTick', [0.1 0.3 1 3 10 30 100]);
                title(ax, sprintf('Rs/Rc: %2.2f, S/CintSens: %2.2f', theCurrentSTFdata.fittedDoGModelRsRcRatio, theCurrentSTFdata.fittedDoGModelSCIntSensRatio))
                drawnow;

                if (~isempty(videoOBJ))
                    videoOBJ.writeVideo(getframe(hFigProgress));
                end
            end
end

function plotProfiles(ax, supportDegs, centerProfile, surroundProfile)

    cMapEntries = 1024;
    cMap = (brewermap(cMapEntries, '*RdBu')).^1.0;

    mxMin = min(supportDegs);
    mxMax = max(supportDegs);

    x1 = ceil(mxMin/0.05)*0.05;
    x2 = floor(mxMax/0.05)*0.05;

    sensitivityRange(1) = -1.05*max(surroundProfile);
    sensitivityRange(2) = 1.05*max(centerProfile);
    sensitivityTicks = -(1:0.1:1)*max(centerProfile);

    shadedAreaBetweenTwoLines(ax, supportDegs, centerProfile, ...
                centerProfile*0, cMap(512+256,:), 'none', 0.3, 1.5, '-');
    hold(ax, 'on');
    shadedAreaBetweenTwoLines(ax, supportDegs, -surroundProfile, ...
                -surroundProfile*0, cMap(512-256,:), 'none', 0.3, 1.5, '-');
    plot(ax,supportDegs, centerProfile-surroundProfile, 'k-', 'LineWidth', 1.0);
    plot(ax,supportDegs, supportDegs*0, 'k--', 'LineWidth', 0.5);
    set(ax, 'XTick', x1: 0.05: x2, 'XLim', [mxMin mxMax], 'YLim', sensitivityRange, 'YTick', sensitivityTicks);
    xlabel(ax, 'space (degs)');
    ylabel(ax, 'gain');
end

function shadedAreaBetweenTwoLines(ax,x,y1, y2, faceColor, edgeColor, faceAlpha, lineWidth, lineStyle)
    x = [x  x(end)  fliplr(x)  x(1)];
    y = [y1 y2(end) fliplr(y2) y2(1)];
    px = reshape(x, [1 numel(x)]);
    py = reshape(y, [1 numel(y)]);
    pz = -10*eps*ones(size(py)); 
    patch(ax,px,py,pz,'FaceColor',faceColor,'EdgeColor', edgeColor, ...
        'FaceAlpha', faceAlpha, 'LineWidth', lineWidth, 'LineStyle', lineStyle);
end
