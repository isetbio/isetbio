function mask = contrastMask(spatialSupport, varargin)
    % Parse optional input
    p = inputParser;
    p.addParameter('zeroInCentralRegionWithSize', [], @(x)((numel(x) == 1)||(numel(x) == 2)));
    p.parse(varargin{:});
    zeroInCentralRegionWithSize = p.Results.zeroInCentralRegionWithSize;

    x = spatialSupport;
    y = spatialSupport;
    [X,Y] = meshgrid(x,y);

    % Default mask: all ones
    mask = X*0+1;

    if (~isempty(zeroInCentralRegionWithSize))
        if (numel(zeroInCentralRegionWithSize) == 1)
            idx = find(abs(X) < 0.5*zeroInCentralRegionWithSize(1));
        else
            idx = find(...
                (abs(X) < 0.5*zeroInCentralRegionWithSize(1)) & ...
                (abs(Y) < 0.5*zeroInCentralRegionWithSize(2)) ...
                );
        end

        mask(idx) = 0;
    end

end

       