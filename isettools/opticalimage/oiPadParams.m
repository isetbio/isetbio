function [padSize, padValue, rectRadius] = oiPadParams(oi)
% Determine padding params from the oi.pad struct
%
% Syntax:
%   [padSize, padValue, rectRadius] = oiPadParams(oi)
%
% Description:
%    Determine the padSize, padValue and rectRadius padding params from the
%    oi.pad struct if the oi.pad has been set, e.g., by the call:
%
%         oi = oiSet(oi, 'pad', padStruct)
%
%    Otherwise, select the defaults:
%       - padSize  = imSize/8
%       - padValue = 'mean photons'
%
% Inputs:
%    oi         - Struct. An optical image structure
%
% Outputs:
%    padSize    - Matrix. A matrix containing the dimensions to pad out.
%    padValue   - String. How to pad. See validatePadStruct() for valid
%                 padValue values.
%    rectRadius - Number. The half-width of the cropped oi in pixels. This
%                 is set to empty if the requested pad size is larger that
%                 the default pad size (1/8 of the oi width)
%
% Optional key/value pairs:
%    None.
%

% History:
%    10/01/18  npc  Wrote it
%    06/24/19  JNM  Documentation pass

    % By default do not crop the oi
    rectRadius = [];

    % Default pad size
    imSize = oiGet(oi, 'size');
    padSizeDefault  = [round(imSize / 8) 0];
    padValueDefault = 'mean photons';

    if ~isfield(oi, 'pad')
        padSize = padSizeDefault;
        padValue = padValueDefault;
        return;
    end

    if isfield(oi.pad, 'sizeDegs')
        % Convert pad size in degs to pad size in pixels
        oiHorizontalFOVdegs = oiGet(oi, 'hfov');
        cols = oiGet(oi, 'cols');
        sampleSizeDegs = oiHorizontalFOVdegs / cols;
        extraCols = ...
            (oi.pad.sizeDegs - oiHorizontalFOVdegs) / sampleSizeDegs;

        oiVerticalFOVdegs = oiGet(oi, 'vfov');
        rows = oiGet(oi, 'rows');
        sampleSizeDegs = oiVerticalFOVdegs / rows;
        extraRows = (oi.pad.sizeDegs - oiVerticalFOVdegs) / sampleSizeDegs;

        % Define padSize (in pixels)
        if (extraCols * 0.5 > padSizeDefault(1)) && ...
                (extraRows * 0.5 > padSizeDefault(2))
            % Larger than default pad size
            padSize = [round(extraRows * 0.5), ...
                round(extraCols * 0.5), padSizeDefault(3)];
        else
            % Smaller pad size than default pad size.
            % Use default but define a cropRadius
            padSize = padSizeDefault;
            rectRadius = round(0.5 * ( oi.pad.sizeDegs / sampleSizeDegs));
        end
    else  % Default behavior
        padSize = padSizeDefault;
    end

    if isfield(oi.pad, 'value')
        padValue = oi.pad.value;
    else  % Default behavior
        padValue = padValueDefault;
    end

end