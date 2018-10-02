function [padSize, padValue] = oiPadParams(oi)
% Determine the padSize, padValue padding params from the oi.pad struct
%
% Syntax:
%   [padSize, padValue] = oiPadParams(oi)
%
% Description:
%    Determine the padSize, padValue padding params from the oi.pad struct
%    if the oi.pad has been set using an oi = oiSet(oi, 'pad', padStruct);
%    Otherwise, select the defaults:
%       - padSize  = imSize/8
%       - padValue = 'meanphotons'
%
% Inputs:
%    oi        - Struct. An optical image structure
%
% Outputs:
%    padSize   - Matrix.  A matrix containing the dimensions to pad out.
%    padValue  - String.  How to pad. See validatePadStruct() for valid padValue values
%
% Optional key/value pairs:
%    None.
%
% History:
%    10/01/18 npc  Wrote it
%

    % Default pad size
    imSize = oiGet(oi, 'size');
    padSizeDefault  = [round(imSize / 8) 0];
    padValueDefault = 'mean photons';
    
    if (~isfield(oi, 'pad'))
        padSize = padSizeDefault;
        padValue = padValueDefault;
        return;
    end
    
    if (isfield(oi.pad, 'sizeDegs'))
        % Convert pad size in degs to pad size in pixels
        oiHorizontalFOVdegs = oiGet(oi, 'hfov');
        cols = oiGet(oi, 'cols');
        sampleSizeDegs = oiHorizontalFOVdegs/cols;
        extraCols = ceil(0.5*(oi.pad.sizeDegs - oiHorizontalFOVdegs)/sampleSizeDegs);
        
        oiVerticalFOVdegs = oiGet(oi, 'vfov');
        rows = oiGet(oi, 'rows');
        sampleSizeDegs = oiVerticalFOVdegs/rows;
        extraRows = ceil(0.5*(oi.pad.sizeDegs - oiVerticalFOVdegs)/sampleSizeDegs);

        
        % pad size in pixels. If new padSize is < pad size, switch to
        % default pad size
        padSize = [...
            max([padSizeDefault(1) extraRows]) ...
            max([padSizeDefault(2) extraCols])...
            padSizeDefault(3) ...
            ];

    else
        % Default behavior
        padSize = padSizeDefault;
    end
    
    if (isfield(oi.pad, 'value'))
        padValue = oi.pad.value;
    else
        % Default behavior
        padValue = padValueDefault;
    end

end