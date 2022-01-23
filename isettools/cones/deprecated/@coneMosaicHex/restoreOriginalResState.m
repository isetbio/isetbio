function restoreOriginalResState(obj)
    % Restore the original Rectangular Grid State
    %
    % Syntax:
    %   restoreOriginalResState(obj)
    %
    % Description:
    %    Restore the original rectangular grid state
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
    obj.patternSampleSize = obj.patternSampleSizeOriginatingRectGrid;
    obj.mosaicSize = size(obj.patternOriginatingRectGrid);
    obj.pattern = obj.patternOriginatingRectGrid;    
    obj.setSizeToFOV(obj.fovOriginatingRectGrid);
end