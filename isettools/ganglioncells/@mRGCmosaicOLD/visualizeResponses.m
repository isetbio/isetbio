function visualizeResponses(obj, responseTemporalSupportSeconds, mRGCMosaicResponse, varargin)
    p = inputParser;
    p.addParameter('stimulusTemporalSupportSeconds', [], @isnumeric);
    p.addParameter('stimulusSceneSequence',[], @(x)(isempty(x) || iscell(x)));
    p.addParameter('opticalSequence',[], @(x)(isempty(x) || isa(x, 'oiSequence') || (isa(x, 'oiArbitrarySequence'))));
    p.addParameter('coneMosaicResponse', []);
    p.addParameter('coneMosaicResponseType', 'excitations', @(x)(ischar(x) && (ismember(x, {'excitations', 'photocurrents'})) ));
    p.addParameter('coneMosaicResponseTemporalSupportSeconds', []);
    p.parse(varargin{:});  
    
    stimulusTemporalSupportSeconds = p.Results.stimulusTemporalSupportSeconds;
    stimulusSequence = p.Results.stimulusSceneSequence;
    opticalSequence = p.Results.opticalSequence;
    coneMosaicResponse = p.Results.coneMosaicResponse;
    coneMosaicResponseType = p.Results.coneMosaicResponseType;
    coneMosaicResponseTemporalSupportSeconds = p.Results.coneMosaicResponseTemporalSupportSeconds;
    
    % Validate inputs
    validateInput(obj, responseTemporalSupportSeconds, mRGCMosaicResponse, ...
        stimulusTemporalSupportSeconds, stimulusSequence, opticalSequence, ...
        coneMosaicResponseTemporalSupportSeconds, coneMosaicResponse);
    
    nInstances = size(mRGCMosaicResponse,1);
    rgcsNum = size(mRGCMosaicResponse,2);
    tBinsNum = size(mRGCMosaicResponse,3);
    
    if (~isempty(coneMosaicResponse))
        % If we get no coneMosaicResponseTemporalSupport, assume it is the
        % same as the mRGC temporal support
        if (isempty(coneMosaicResponseTemporalSupportSeconds))
            coneMosaicResponseTemporalSupportSeconds = responseTemporalSupportSeconds;
        end
        conesNum = size(coneMosaicResponse,2);
    end
    
    meanConeMosaicResponse = squeeze(mean(coneMosaicResponse,1));
    meanMRGCResponse = squeeze(mean(mRGCMosaicResponse,1));
    
    hFig = figure(100);
    if (isempty(coneMosaicResponse))
        set(hFig, 'Position', [10 10 800 900]);
        stimAxes = subplot(2,2,1);
        stimLuminanceTracesAxes = subplot(2,2,2);
        mRGCMosaicAxes = subplot(2,2,3);
        mRGCTraceAxes = subplot(2,2,4);
        coneMosaicAxes = [];
        coneTraceAxes = [];
    else
        set(hFig, 'Position', [10 10 800 1300]);
        stimAxes = subplot(3,2,1);
        stimLuminanceTracesAxes = subplot(3,2,2);
        coneMosaicAxes = subplot(3,2,3);
        coneTraceAxes = subplot(3,2,4);
        mRGCMosaicAxes = subplot(3,2,5);
        mRGCTraceAxes = subplot(3,2,6);
    end
    
    stimLuminanceTrace = zeros(1,numel(stimulusTemporalSupportSeconds));
    coneTrace = zeros(1, numel(responseTemporalSupportSeconds));
    mRGCTrace = zeros(1, numel(responseTemporalSupportSeconds));
    

    stimFrameInterval = stimulusTemporalSupportSeconds(2)-stimulusTemporalSupportSeconds(1);
    coneResponseRange = prctile(meanConeMosaicResponse(:), [5 95]);
    mRGCResponseRange = [min(meanMRGCResponse(:)) max(meanMRGCResponse(:))];
    
    for tBin = 1:numel(responseTemporalSupportSeconds)
        time = responseTemporalSupportSeconds(tBin);
        stimulusTimeBin = floor(time/stimFrameInterval)+1;
        
        % Update the current stimulus
        stimLuminanceTrace(tBin) = renderCurrestStimulus(stimAxes, stimulusSequence{stimulusTimeBin}, time);
        
        % Update the stimulus luminance trace plot
        plot(stimLuminanceTracesAxes, responseTemporalSupportSeconds(1:tBin)*1000, stimLuminanceTrace(1:tBin), 'ks-');
        set(stimLuminanceTracesAxes, 'XTick', 0:100:1000, 'FontSize', 18);
        axis(stimLuminanceTracesAxes, 'square');
        title(stimLuminanceTracesAxes,'Stimulus luminance')
        xlabel(stimLuminanceTracesAxes, 'time (msec)');
        ylabel(stimLuminanceTracesAxes,'Cd/m2');
        
        if (~isempty(coneMosaicAxes))
            % Update the cone mosaic response plot
            coneTrace(tBin) = renderCurrentConeMosaicActivation(coneMosaicAxes, obj, squeeze(meanConeMosaicResponse(:, tBin)), ...
                coneResponseRange, coneMosaicResponseType);
            % Update the cone response trace plot
            plot(coneTraceAxes, responseTemporalSupportSeconds(1:tBin)*1000, coneTrace(1:tBin), 'ks-');
            set(coneTraceAxes, 'XTick', 0:100:1000, 'FontSize', 18);
            title(coneTraceAxes, 'Single cone response');
            xlabel(coneTraceAxes, 'time (msec)');
            switch (coneMosaicResponseType)
                case 'excitations'
                    ylabel(coneTraceAxes, 'R*/sec');
                case 'photocurrents'
                    ylabel(coneTraceAxes, 'pAmps');
            end
        end
        
        % Update the mRGC mosaic plot
        mRGCTrace(tBin) = renderCurrentMRGCMosaicActivation(mRGCMosaicAxes, obj, squeeze(meanMRGCResponse(:, tBin)), mRGCResponseRange);

        % Update the mRGC response trace plot
        plot(mRGCTraceAxes, responseTemporalSupportSeconds(1:tBin)*1000, mRGCTrace(1:tBin), 'ks-');
        set(mRGCTraceAxes , 'XTick', 0:100:1000, 'FontSize', 18);
        axis(mRGCTraceAxes, 'square');
        title(mRGCTraceAxes, 'Single mRGC response');
        xlabel(mRGCTraceAxes, 'time (msec)');
        drawnow;
    end
    
    
    hFig = figure(1000);
    set(hFig, 'Name', 'spatiotemporal cone and mRGC mosaic responses');
    if (~isempty(coneMosaicResponse))
        % Plot the mean cone spatiotemporal response
        subplot(2,1,1);
        imagesc(coneMosaicResponseTemporalSupportSeconds, 1:conesNum, meanConeMosaicResponse);
        title('cone mosaic response');
        xlabel('time (seconds)');
        ylabel('cone number');
    end
     
    % Plot the mean mRGC spatiotemporal response
    subplot(2,1,2);
    imagesc(responseTemporalSupportSeconds, 1:rgcsNum, meanMRGCResponse);
    title('mRGC mosaic response');
    xlabel('time (seconds)');
    ylabel('mRGC cell number');
    colormap(gray(1024));
end

function luminanceTrace = renderCurrestStimulus(ax, stimulusScene, time)
    cla(ax);
    lumImage = sceneGet(stimulusScene, 'luminance');
    m = floor(size(lumImage,1)/2);
    n = floor(size(lumImage,2)/2);
    luminanceTrace = lumImage(m,n);
    stimSize = sceneGet(stimulusScene, 'size');
    xFOVdegs = sceneGet(stimulusScene, 'wAngular');
    yFOVdegs = sceneGet(stimulusScene, 'hAngular');
    xDegs = linspace(0,xFOVdegs,stimSize(2));
    yDegs = linspace(0,yFOVdegs,stimSize(1));
    xDegs = xDegs - mean(xDegs);
    yDegs = yDegs - mean(yDegs);
    imagesc(ax, xDegs, yDegs, sceneGet(stimulusScene, 'rgbimage'));
    xlabel(ax, 'space (degs)');
    ylabel(ax, 'space (degs)');
    axis(ax, 'image');
    set(ax, 'XTick', -10:0.25:10, 'YTick', -10:0.25:10, 'FontSize', 18);
    title(ax, sprintf('scene @ %2.0f msec', time*1000));
end

function coneTrace = renderCurrentConeMosaicActivation(ax, mRGCMosaicOBJ, coneMosaicResponse, signalRange, signalType)
    cla(ax)
    tracedPositionDegs = [0 0];
    cmStruct = mRGCMosaicOBJ.inputConeMosaic.geometryStructAlignedWithSerializedConeMosaicResponse();
    [~, tracedConeIndex] = min(sum((bsxfun(@minus, cmStruct.coneLocs, tracedPositionDegs)).^2,2));
    switch signalType
        case 'excitations'
            coneTrace = coneMosaicResponse(tracedConeIndex,:) / mRGCMosaicOBJ.inputConeMosaic.integrationTime;
        case 'photocurrents'
            coneTrace = coneMosaicResponse(tracedConeIndex,:);
    end
    
    mRGCMosaicOBJ.inputConeMosaic.renderActivationMap(ax, coneMosaicResponse, ...
        'signalRange', signalRange, ...
        'visualizedFOV', 1.0, ...
        'mapType', 'modulated disks', ...
	    'backgroundColor', [0 0 0]);
    title(ax, 'Cone mosaic activation');
end

function mRGCTrace = renderCurrentMRGCMosaicActivation(ax, mRGCMosaicOBJ, mRGCMosaicResponse, signalRange)
    
    tracedPositionDegs = [0 0];
    center = mean(mRGCMosaicOBJ.rgcRFpositionsDegs,1);
    targetPos = tracedPositionDegs+center;
    diffs = bsxfun(@minus, mRGCMosaicOBJ.rgcRFpositionsDegs, targetPos);
    distancesToTracePosition = sqrt(sum(diffs.^2,2));
    [~,idx] = min(distancesToTracePosition);
    mRGCTrace = mRGCMosaicResponse(idx);
    
    cla(ax)
    mRGCMosaicOBJ.visualizeActivationMap(ax, mRGCMosaicResponse, ...
        'signalRange', signalRange, ...
        'domain', 'degs', ...
        'visualizedFOV', [1 1]);
    title(ax, 'mRGC mosaic activation');

end


function validateInput(obj, responseTemporalSupportSeconds, mRGCMosaicResponse, ...
    stimulusTemporalSupportSeconds, stimulusSequence, opticalSequence, ...
    coneMosaicResponseTemporalSupportSeconds, coneMosaicResponse)

    % Assert that the dimensionality of mRGCresponse is consistent with the temporal support 
    assert(size(responseTemporalSupportSeconds,2) == size(mRGCMosaicResponse,3), ...
        'mismatch in timeBins between responseTemporalSupport and mRGCMosaicResponse.');

    % Assert that the #of RGCs in the mRGC mosaic is constistent with
    % the mRGCMosaicResponse
    assert(size(obj.rgcRFpositionsMicrons,1) == size(mRGCMosaicResponse,2), ...
            'mismatch in number of RGCs in mRGCMosaic and mRGCMosaicResponse');
        
    if (~isempty(stimulusSequence)) 
        % Assert that the stimulusTemporalSupport vector is non-empty
        assert(~isempty(stimulusTemporalSupportSeconds), ...
            '''stimulusSequence'' was passed without ''stimulusTemporalSupport'' information.');
        % Assert that the dimensionality of stimulusSequence matches that of stimulusTemporalSupportSeconds
        assert(numel(stimulusSequence) == numel(stimulusTemporalSupportSeconds), ...
            'mismatch in stimulus sequence and stimulusTemporalSupport dimensionalities');
    end
    
    if (~isempty(opticalSequence)) && (~isempty(stimulusSequence)) 
        % Assert that the dimensionality of the optical sequence matches
        % that of the stimulus sequence
        assert(numel(opticalSequence.timeAxis) == numel(stimulusSequence), ...
            'mismatch in stimulus and optical sequence dimensionalities');
    end
        
    if (~isempty(coneMosaicResponse))
        % If we get no coneMosaicResponseTemporalSupport, assume it is the
        % same as the mRGC temporal support
        if (isempty(coneMosaicResponseTemporalSupportSeconds))
            coneMosaicResponseTemporalSupportSeconds = responseTemporalSupportSeconds;
        end
        
        % Assert that the dimensionality of coneMosaicResponse is consistent with the temporal support 
        assert(size(coneMosaicResponseTemporalSupportSeconds,2) == size(coneMosaicResponse,3), ...
            'mismatch in timeBins between responseTemporalSupport and coneMosaicResponse.');
    
        % Assert that the #of response instances match
        assert(size(coneMosaicResponse,1) == size(mRGCMosaicResponse,1), ...
            'mismatch in number of instances between mRGC mosaic responses and cone mosaic responses');
        
        % Assert that the #of cones in the cone mosaic is constitent with
        % the coneMosaicResponse
        assert(size(obj.inputConeMosaicMetaData.conePositionsDegs,1) == size(coneMosaicResponse,2), ...
            'mismatch in number of cones in coneMosaic and coneMosaicResponse');
        
    end
    
end

