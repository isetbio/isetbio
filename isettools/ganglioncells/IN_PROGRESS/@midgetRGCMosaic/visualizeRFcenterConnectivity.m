function visualizeRFcenterConnectivity(obj, varargin)
    p = inputParser;
    p.addParameter('figureHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('axesHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('titleString', '', @(x)(isempty(x) || (ischar(x))));
    p.parse(varargin{:});
    
    hFig = p.Results.figureHandle;
    ax = p.Results.axesHandle;
    titleString = p.Results.titleString;

    obj.theMosaicConnectorOBJ.visualizeDestinationLatticePooling(...
        'figureHandle', hFig, ...
        'axesHandle', ax, ...
        'titleString', titleString, ...
        'titleWithPoolingStats', true);
end