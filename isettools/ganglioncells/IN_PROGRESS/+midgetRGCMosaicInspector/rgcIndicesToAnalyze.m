function rgcIndices = rgcIndicesToAnalyze(mosaicFileName, varargin)
    % Parse input
    p = inputParser;
    p.addParameter('maxRGCsNum', [], @(x)(isempty(x) || (isscalar(x))));
    p.parse(varargin{:});

    maxRGCsNum = p.Results.maxRGCsNum;


    % Load the fully-connected mosaic
    load(mosaicFileName, 'theMidgetRGCmosaic');

    rgcsNum = size(theMidgetRGCmosaic.rgcRFsurroundConePoolingMatrix,2);
    if (isempty(maxRGCsNum))
        skippedRGCs = 1;
    else
        skippedRGCs = max([1 round(rgcsNum/maxRGCsNum)]);
    end

    % Sort RGCs with respect to their distance from the mRGC mosaic center
    mRGCmosaicCenterDegs = mean(theMidgetRGCmosaic.rgcRFpositionsDegs,1);
    d = sum((bsxfun(@minus, theMidgetRGCmosaic.rgcRFpositionsDegs, mRGCmosaicCenterDegs)).^2,2);
    [~,sortedRGCindices] = sort(d, 'ascend');

    rgcIndices = sortedRGCindices(1:skippedRGCs:rgcsNum);
end