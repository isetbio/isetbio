function saveOriginalResState(obj)
    % Save the original Rectangular Grid State
    %
    % Syntax:
    %   saveOriginalResState(obj, fov)
    %
    % Description:
    %    Save the original rectangular grid state. Internal method used
    %    only by the setSizeToFOVForHexMosaic.
    %
    % Inputs:
    %    obj - The cone mosaic hex object
    %
    % Outputs:
    %    None.
    %
    % Optional key/value pairs:
    %    None.
    %
    obj.fovOriginatingRectGrid = obj.fov;
    obj.coneLocsOriginatingRectGrid = obj.coneLocs;
    obj.patternOriginatingRectGrid = obj.pattern;
    obj.patternSampleSizeOriginatingRectGrid = obj.patternSampleSize;
end
