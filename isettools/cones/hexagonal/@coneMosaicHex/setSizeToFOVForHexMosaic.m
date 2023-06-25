function setSizeToFOVForHexMosaic(obj, fov)
    % Set the size to the field of view for the hex mosaic
    %
    % Syntax:
    %   setSizeToFOVForHexMosaid(obj, fov)
    %
    % Description:
    %    Set the size to the field of view (FOV) for the hex mosaic.
    %
    % Inputs:
    %    obj - The cone mosaic hex object
    %    fov - The field of view
    %
    % Outputs:
    %    None.
    %
    % Optional key/value pairs:
    %    None.
    %
    obj.restoreOriginalResState();
    obj.setSizeToFOV(fov);
    obj.saveOriginalResState();
    obj.resampleGrid(obj.resamplingFactor);
end
