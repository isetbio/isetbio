function resampleGrid(obj, resamplingFactor)
%  % Sample the original rectangular mosaic using a hex grid sampled at the passed resamplingFactor
%
% NPC, ISETBIO TEAM, 2015

    % Restore original state
    obj.restoreOriginalResState();
    
    % Compute perfect hex grid
    obj.coneLocsHexGrid = computeHexGridNodes(obj);
    
    % Decrease patternSampleSize and increase mosaicSize both by the resamplingFactor
    obj.resamplingFactor = resamplingFactor;
    obj.patternSampleSize = obj.patternSampleSize / obj.resamplingFactor;
    obj.mosaicSize = obj.mosaicSize * obj.resamplingFactor;
    
    % Sample the perfect hex grid on the high res rectangular grid nodes
    obj.pattern = rectSampledHexPattern(obj);
    
    % Plot grid
    obj.visualizeGrid('generateNewFigure', true, 'panelPosition', [1 1]);
end

function hexLocs = computeHexGridNodes(obj)

    fprintf('Computing hex mosaic corresponding to rect mosaic with %d rows and %d cols\n', obj.rows, obj.cols);
    extraCols = round(obj.cols/(sqrt(3)/2)) - obj.cols;
    rectXaxis2 = (1:(obj.cols+extraCols));
    [X2,Y2] = meshgrid(rectXaxis2, 1:obj.rows);
    
    X2 = X2 * sqrt(3)/2;
    for iCol = 1:size(Y2,2)
        Y2(:,iCol) = Y2(:,iCol) - mod(iCol-1,2)*0.5;
    end
    
    densityScalar = sqrt(2/sqrt(3.0));
    X2 = X2 * densityScalar;
    Y2 = Y2 * densityScalar;
    marginInConePositions = 0.1;
    indicesToKeep = (X2>=-marginInConePositions) & ...
                    (X2<=obj.cols+marginInConePositions) &...
                    (Y2>=-marginInConePositions) & ...
                    (Y2<=obj.rows+marginInConePositions);             
    xHex = X2(indicesToKeep);
    yHex = Y2(indicesToKeep);
   
    xHex = xHex * obj.patternSampleSize(1);
    yHex = yHex * obj.patternSampleSize(2);
    hexLocs = [xHex(:)-mean(xHex(:))   yHex(:)-mean(yHex(:))];
end


function pattern = rectSampledHexPattern(obj)

    xHex = squeeze(obj.coneLocsHexGrid(:,1));
    yHex = squeeze(obj.coneLocsHexGrid(:,2));
    
    xRectHiRes = (1:obj.cols) * obj.patternSampleSize(1); xRectHiRes = xRectHiRes - mean(xRectHiRes);
    yRectHiRes = (1:obj.rows) * obj.patternSampleSize(2); yRectHiRes = yRectHiRes - mean(yRectHiRes);
    [xx,yy] = meshgrid(xRectHiRes, yRectHiRes); 
    
    % Optimize sub-sampling jitter. In progress...
    optimizeSubSamplingJitter = false;
    if (optimizeSubSamplingJitter)
        [xHex, yHex] = optimizeJitter(xHex, yHex, xx, yy);
    end
    
    xRectOriginal = (1:size(obj.patternOriginatingRectGrid,2)) * obj.patternSampleSizeOriginatingRectGrid(1); xRectOriginal = xRectOriginal - mean(xRectOriginal);
    yRectOriginal = (1:size(obj.patternOriginatingRectGrid,1)) * obj.patternSampleSizeOriginatingRectGrid(2); yRectOriginal = yRectOriginal - mean(yRectOriginal);
    xRectOriginal = xRectOriginal / max(xRectOriginal) * max(xRectHiRes);
    yRectOriginal = yRectOriginal / max(yRectOriginal) * max(yRectHiRes);
    [xxx,yyy] = meshgrid(xRectOriginal, yRectOriginal);
    
    pattern = zeros(numel(yRectHiRes), numel(xRectHiRes))+1;    
    for k = 1:numel(xHex)
        % determine cone ID of closest cone in original pattern
        [~, idx] = min(sqrt((xHex(k)-xxx(:)).^2 + (yHex(k)-yyy(:)).^2));
        [originalRow, originalCol] = ind2sub(size(obj.patternOriginatingRectGrid), idx);
        coneID = obj.patternOriginatingRectGrid(originalRow,originalCol);
        
        [~, idx] = min(sqrt((xHex(k)-xx(:)).^2 + (yHex(k)-yy(:)).^2));
        [hiResPatternRow,hiResPatternCol] = ind2sub(size(pattern), idx);
        pattern(hiResPatternRow, hiResPatternCol) = coneID;
    end
    
    debugPlots = false;
    if (debugPlots)
        cmap = [0 0 0; 1 0 0; 0 1 0; 0 0 1];
        figure(10); clf;
        imagesc((obj.patternOriginatingRectGrid-1)/3);
        axis 'equal'; axis 'xy'
        set(gca, 'CLim', [0 1]);
        colormap(cmap);
        axis 'image'

        figure(11); clf;
        imagesc((pattern-1)/3);
        axis 'equal'; axis 'xy'
        set(gca, 'CLim', [0 1]);
        colormap(cmap);
        axis 'image'
    end
end

function [xHex, yHex] = optimizeJitter(xHex, yHex, xx, yy)
    tic
    fprintf('\n Optimizing resampling grid jitter ...\n');
    
    scalings = 1 + linspace(0.0,0.05,15);
    jitterListX = scalings;
    jitterListY = jitterListX;
    
    gridPointIndicesToCheck = find(xHex >= -Inf & yHex >= -Inf); 
    xHexTmp = xHex(gridPointIndicesToCheck);
    yHexTmp = yHex(gridPointIndicesToCheck);
    
    errors = zeros(numel(jitterListY), numel(jitterListX), numel(xHexTmp));
    for yk = 1:numel(jitterListY)
        yTmp = yHexTmp * jitterListY(yk);
        for xk = 1:numel(jitterListX)
            xTmp = xHexTmp * jitterListX(xk);
            for k = 1:numel(xHexTmp)
                errors(yk,xk,k) = min(sqrt((xTmp(k)-xx(:)).^2 + (yTmp(k)-yy(:)).^2));
            end
        end
    end
    
    % max across all points
    totalError = sum(errors, 3);
    
    % find jitter coords that result in min (totalError)
    [~, idx] = min(totalError(:));
    [ykBest, xkBest] = ind2sub(size(totalError), idx);
    fprintf('Done in %2.1f seconds\n', toc)
    
    figure(111); clf;
    imagesc(jitterListX, jitterListY, totalError);
    axis 'equal'; axis 'xy';
    hold on;
    plot(jitterListX(xkBest), jitterListY(ykBest), 'ws');
    drawnow;
    
    xHex = xHex * jitterListX(xkBest);
    yHex = yHex * jitterListY(ykBest);
end
