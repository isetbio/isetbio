function [fixationMap, fixationMapSupportX, fixationMapSupportY, ...
        fixationMapXSlice, fixationMapYSlice] = computeFixationMap(timeAxis, emPaths, emPosRange, emPosDelta, varargin)
    
    p = inputParser;
    addParameter(p, 'maxDurationSeconds', Inf, @isnumeric);
    parse(p, varargin{:});
    
    % Only analyze span within the maxDurationSeconds
    idx = find(timeAxis <= p.Results.maxDurationSeconds);
    emPaths = emPaths(:,idx,:);
    
    xEdges = emPosRange(1):emPosDelta:emPosRange(end)+emPosDelta;
    yEdges = emPosRange(1):emPosDelta:emPosRange(end)+emPosDelta;
    xPos = squeeze(emPaths(:,:,1)); yPos = squeeze(emPaths(:,:,2));
    fixationMap = histcounts2(xPos(:),yPos(:),xEdges, yEdges);
    fixationMap = fixationMap/max(fixationMap(:));
    fixationMapSupportX = xEdges(1:end-1)+emPosDelta/2;
    fixationMapSupportY = yEdges(1:end-1)+emPosDelta/2;
    
    [~,midRow] = min(abs(yEdges));
    fixationMapXSlice = squeeze(fixationMap(midRow,:));
    fixationMapYSlice = squeeze(fixationMap(:,midRow));        
end