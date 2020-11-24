function visualizeResponses(obj, responseTemporalSupportSeconds, mRGCMosaicResponse, varargin)
    p = inputParser;
    p.addParameter('stimulusTemporalSupportSeconds', [], @isnumeric);
    p.addParameter('stimulusSceneSequence',[], @(x)(isempty(x) || iscell(x)));
    p.addParameter('opticalSequence',[], @(x)(isempty(x) || isa(x, 'oiSequence') || (isa(x, 'oiArbitrarySequence'))));
    p.addParameter('coneMosaicResponse', []);
    p.addParameter('coneMosaicResponseTemporalSupportSeconds', []);
    p.parse(varargin{:});  
    
    stimulusTemporalSupportSeconds = p.Results.stimulusTemporalSupportSeconds;
    stimulusSequence = p.Results.stimulusSceneSequence;
    opticalSequence = p.Results.opticalSequence;
    coneMosaicResponse = p.Results.coneMosaicResponse;
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
    set(hFig, 'Position', [10 10 800 1300]);
    stimAxes = subplot(3,2,1);
    coneMosaicAxes = subplot(3,2,3);
    mRGCMosaicAxes = subplot(3,2,5);
    stimTracesAxes = subplot(3,2,2);
    coneTraceAxes = subplot(3,2,4);
    mRGCTraceAxes = subplot(3,2,6);
    
    stimTrace = zeros(1,numel(stimulusTemporalSupportSeconds));
    coneTrace = zeros(1, numel(responseTemporalSupportSeconds));
    mRGCTrace = zeros(1, numel(responseTemporalSupportSeconds));
    
    tracedPositionDegs = [0 0];
    
    stimFrameInterval = stimulusTemporalSupportSeconds(2)-stimulusTemporalSupportSeconds(1);
    coneResponseRange = prctile(meanConeMosaicResponse(:), [5 95]);
    mRGCResponseRange = [min(meanMRGCResponse(:)) max(meanMRGCResponse(:))];
    
    for tBin = 1:numel(responseTemporalSupportSeconds)
        time = responseTemporalSupportSeconds(tBin);
        stimulusTimeBin = floor(time/stimFrameInterval)+1;
        stimTrace(tBin) = renderCurrestStimulus(stimAxes, stimulusSequence{stimulusTimeBin}, time, tracedPositionDegs);
        coneTrace(tBin) = renderCurrentConeMosaicActivation(coneMosaicAxes, obj, squeeze(meanConeMosaicResponse(:, tBin)), coneResponseRange, tracedPositionDegs);
        mRGCTrace(tBin) = renderCurrentMRGCMosaicActivation(mRGCMosaicAxes, obj, squeeze(meanMRGCResponse(:, tBin)), mRGCResponseRange, tracedPositionDegs);
        subplot(3,2,2);
        plot(responseTemporalSupportSeconds(1:tBin)*1000, stimTrace(1:tBin), 'ks-');
        set(gca, 'XTick', 0:50:1000);
        subplot(3,2,4);
        plot(responseTemporalSupportSeconds(1:tBin)*1000, coneTrace(1:tBin), 'ks-');
         set(gca, 'XTick', 0:50:1000);
        subplot(3,2,6);
        plot(responseTemporalSupportSeconds(1:tBin)*1000, mRGCTrace(1:tBin), 'ks-');
         set(gca, 'XTick', 0:50:1000);
        drawnow;
    end
    
    
    figure(1000);
    
    if (~isempty(coneMosaicResponse))
        % Plot the mean cone spatiotemporal response
        subplot(2,1,1);
        imagesc(coneMosaicResponseTemporalSupportSeconds, 1:conesNum, meanConeMosaicResponse);
        title('cone mosaic response');
        xlabel('time (seconds)');
        ylabel('cone index');
    end
     
    % Plot the mean mRGC spatiotemporal response
    subplot(2,1,2);
    imagesc(responseTemporalSupportSeconds, 1:rgcsNum, meanMRGCResponse);
    title('mRGC mosaic response');
    xlabel('time (seconds)');
    ylabel('mRGC index');
end

function stimTrace = renderCurrestStimulus(ax, stimulusScene, time, tracedPositionDegs)
    cla(ax);
    lumImage = sceneGet(stimulusScene, 'luminance');
    m = floor(size(lumImage,1)/2);
    n = floor(size(lumImage,2)/2);
    stimTrace = lumImage(m,n);
    imagesc(ax, sceneGet(stimulusScene, 'rgbimage'));
    axis(ax, 'image');
    set(ax, 'XTick', [], 'YTick', []);
    title(ax, sprintf('%d msec', time*1000));
end

function coneTrace = renderCurrentConeMosaicActivation(ax, mRGCMosaicOBJ, coneMosaicResponse, signalRange, tracedPositionDegs)
    cla(ax)
    coneTrace = 0;
    mRGCMosaicOBJ.inputConeMosaic.renderActivationMap(ax, coneMosaicResponse, ...
        'signalRange', signalRange, ...
        'visualizedFOV', 1.0, ...
        'mapType', 'modulated disks', ...
	    'backgroundColor', [0 0 0]);
    
end

function mRGCTrace = renderCurrentMRGCMosaicActivation(ax, mRGCMosaicOBJ, mRGCMosaicResponse, signalRange, tracedPositionDegs)
    
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

