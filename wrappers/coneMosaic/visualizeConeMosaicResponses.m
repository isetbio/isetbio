function visualizeConeMosaicResponses(coneMosaic, responses, responseSignalName, varargin)
p = inputParser;
p.addParameter('customTitle', '', @ischar);
% Parse input
p.parse(varargin{:});
customTitle = p.Results.customTitle;

    if (ndims(responses) == 4)
        % Compute the mean differential response
        meanResponse = squeeze(mean(responses, 1));
        meanDiffResponse = bsxfun(@minus, meanResponse, squeeze(meanResponse(:,:,1)));
        % Find peak cone position and peak response time based on mean
        % differential reseponse
        [~,idx] = max(meanDiffResponse(:));
        [rowConeOfPeakResponse,colConeOfPeakResponse,timePointOfPeakResponse] = ...
            ind2sub(size(meanResponse), idx);
        
        patternSupport = coneMosaic.patternSupport;
        peakConePositionMicrons = [...
            patternSupport(rowConeOfPeakResponse,colConeOfPeakResponse,1) ...
            patternSupport(rowConeOfPeakResponse,colConeOfPeakResponse,2)] * 1e6;
        
        singleConeTemporalResponse = squeeze(...
            responses(:,rowConeOfPeakResponse,...
            colConeOfPeakResponse,:));
        
        timeAxis = coneMosaic.timeAxis;
        visualizePeakConeResponseInTime(timeAxis, singleConeTemporalResponse, responseSignalName, peakConePositionMicrons);
        
        % Visualize cone mosaic response at the peak time
        responses = squeeze(responses(:,:,:,timePointOfPeakResponse));
        time = timeAxis(timePointOfPeakResponse);
        responseRange = [min(meanResponse(:)) max(meanResponse(:))];
    else
        time = 0;
        responseRange = [min(responses(:)) max(responses(:))];
    end
    
   
    visualize2DResponseAtASingleTimePoint(coneMosaic, responses, responseRange, responseSignalName, time, customTitle);
end

function visualizePeakConeResponseInTime(timeAxis, peakConeTemporalResponses, signalName, peakConePositionMicrons)
    figure(); clf;
    plot(timeAxis, peakConeTemporalResponses, 'c-', 'LineWidth', 2.0);
    hold on;
    plot(timeAxis, mean(peakConeTemporalResponses,1), 'b-', 'LineWidth', 5.0);
    xlabel('\it time (seconds)');
    ylabel(sprintf('\\it %s', signalName));
    title(sprintf('Cone at %2.2f, %2.2f microns', peakConePositionMicrons(1),peakConePositionMicrons(2)));
    set(gca, 'FontSize', 18);
    grid on; box on
end

function visualize2DResponseAtASingleTimePoint(coneMosaic, responses, responseRange, responseSignalName, responseTime, customTitle)
    figure(); clf;
    subplotRows = 2;
    subplotCols = 3;
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', subplotRows, ...
       'colsNum', subplotCols, ...
       'heightMargin',  0.01, ...
       'widthMargin',    0.02, ...
       'leftMargin',     0.02, ...
       'rightMargin',    0.12, ...
       'bottomMargin',   0.04, ...
       'topMargin',      0.01);

    instancesNum = size(responses,1);
    visualizedTrialsNum = min([subplotRows*subplotCols, instancesNum]);
    
    xTicksDegs =  round([-1 -0.5 0 0.5 1]*coneMosaic.fov(2)/2*100)/100;
    xTicksMeters = xTicksDegs*coneMosaic.micronsPerDegree*1e-6;
    xTickLabelDegs = sprintf('%2.1f\n', xTicksDegs);
    
    axHandle = subplot('Position', subplotPosVectors(1,1).v);
    coneMosaic.visualizeGrid(...
        'axesHandle', axHandle, ...
        'ticksInVisualDegs', true);
    axis(axHandle, 'square')
    set(axHandle, 'FontSize', 12, 'YTickLabel', {}, ...
        'XTick', xTicksMeters, 'XTickLabel', xTickLabelDegs);
    xlabel(axHandle, 'space (degs)')
    ylabel('');
    
    for k = 1:2
        axHandle = subplot('Position', subplotPosVectors(1,k+1).v);
        if (k == 2)
            coneMosaic.renderActivationMap(axHandle, squeeze(responses(k,:,:)), ...
                'mapType', 'modulated disks', ...
                'signalRange', responseRange, ...
                'showColorBar', true, ...
                'labelColorBarTicks', true, ...
                'titleForColorBar', responseSignalName);
        else 
            coneMosaic.renderActivationMap(axHandle, squeeze(responses(k,:,:)), ...
                'mapType', 'modulated disks', ...
                'signalRange', responseRange, ...
                'showColorBar', ~true, ...
                'labelColorBarTicks', ~true);
        end
        ylabel(''); xlabel('');
        set(gca, 'XTick', [], 'YTick', []);
        set(gca, 'FontSize', 12);
        title(sprintf('trial #%d (t: %2.2f sec)', k, responseTime));
    end
    
    % Retrieve indices of cones along horizontal meridian
    [indicesOfConesAlongXaxis, xCoordsOfConesAlongXaxis, typesOfConesAlongXaxis] =...
        indicesOfConesAlongHorizontalMeridian(coneMosaic);

    % Extract the excitations of cones along horizontal meridian
    responses = reshape(responses, [instancesNum  size(responses,2)* size(responses,3)]);
    responsesNxXY = zeros(instancesNum, numel(indicesOfConesAlongXaxis));
    for k = 1:numel(indicesOfConesAlongXaxis)
        responsesNxXY(:,k) = responses(:,indicesOfConesAlongXaxis(k));
    end
    
    % Plot the excitations separately for L-,M- and S-cones
    subplot('Position', [0.08 0.1 0.9 0.33]);
    idx = find(typesOfConesAlongXaxis == 2);
    LconesNum = numel(idx);
    plot(xCoordsOfConesAlongXaxis(idx), responsesNxXY(:,idx), 'r.');
    hold on;
    idx = find(typesOfConesAlongXaxis == 3);
    MconesNum = numel(idx);
    plot(xCoordsOfConesAlongXaxis(idx), responsesNxXY(:,idx), 'g.');
    idx = find(typesOfConesAlongXaxis == 4);
    SconesNum = numel(idx);
    plot(xCoordsOfConesAlongXaxis(idx), responsesNxXY(:,idx), 'b.');
    grid on
    set(gca, 'FontSize', 16, 'XTick', xTicksDegs);
    set(gca, 'YLim', responseRange);
    xlabel('\it space (degs)');
    ylabel(sprintf('\\it %s', responseSignalName));
    if (isempty(customTitle))
        title(gca, sprintf('%d trials, responses of %d L- %d M- and %d-S cones (horiz. meridian)', ...
       	instancesNum, LconesNum, MconesNum, SconesNum), 'FontWeight', 'Normal', 'FontSize', 10);
    else
        title(gca, sprintf('%s - %d trials (cones along horiz. meridian)',customTitle, instancesNum), 'FontWeight', 'Normal', 'FontSize', 10);
    end
end

function [indicesOfConesAlongXaxis,coneXcoordsAlongXaxis, theConeTypes] = indicesOfConesAlongHorizontalMeridian(theMosaic)
    
    conesXcoords = squeeze(theMosaic.patternSupport(1, :, 1));
    conesYcoords = squeeze(theMosaic.patternSupport(:, 1, 2));
    [coneXcoordsMesh, coneYcoordsMesh] = meshgrid(conesXcoords, conesYcoords);
    coneXcoordsDegsMesh = coneXcoordsMesh * 1e6 / theMosaic.micronsPerDegree;
    coneYcoordsDegsMesh = coneYcoordsMesh * 1e6 / theMosaic.micronsPerDegree;
    idx = find(theMosaic.pattern > 1);
    coneXcoordsDegsMesh = coneXcoordsDegsMesh(idx);
    coneYcoordsDegsMesh = coneYcoordsDegsMesh(idx);
    dx = diameterForCircularApertureFromWidthForSquareAperture(...
            theMosaic.pigment.width) * 1e6 / theMosaic.micronsPerDegree;
    indicesOfConesAlongXaxis = find(abs(coneYcoordsDegsMesh) < dx);
    coneXcoordsAlongXaxis = coneXcoordsDegsMesh(indicesOfConesAlongXaxis);
    theConeTypes = theMosaic.pattern(idx(indicesOfConesAlongXaxis));
    indicesOfConesAlongXaxis = idx(indicesOfConesAlongXaxis);
end
