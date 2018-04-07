function DAC = ieLUTLinear(RGB, invGammaTable)
% Convert linear RGB values through a gamma table to DAC values
%
% Syntax:
%   DAC = ieLUTLinear(RGB, invGammaTable)
%
% Description:
%    Convert linear RGB values through a gamma table to DAC values.  This
%    operation is often called gamma correction, and is used to produce the
%    values we pass to display hardware so as to produce desired linear RGB
%    values.
%
%    The input RGB values are assumed to be in the range of [0, 1].
%    They are assumed to be linear with respect to radiance (intensity).
%
%    The input RGB should be an M x N x nPrimaries 3D matrix. However, if
%    the input is NxnPrimaries or it is an nPrimaries vector, we handle
%    that case.
%
%    The size of returned DAC will always be the same as input RGB.
%
%    The returned DAC values are digital values with a bit depth that is
%    determined by the entries in the invGammaTable. 
%
%    We define
%     * The gamma table maps the digital values to the display intensity.
%     * The inverse gamma table maps the display intensity to the digital
%       values.
%
%     The gamma table is directly stored in display calibration files. And
%     the inverse gamma table can be computed via ieLUTInvert.
%
%    If the invGammaTable is a single number, we raise the data to the
%    power invGammaTable.
%
% Inputs:
%    RGB           - Linear RGB Values. 0 - 1. Linear with respect to
%                    radiance (intensity). 
%    invGammaTable - (Optional) Inverse gamma table. This maps the display
%                    intensity to the digital values. If a single number,
%                    raise data to that power. Default is 2.2.
%
% Outputs:
%    DAC           - Converted DAC Values. Follows input RGB's sizing.
%
% Optional key/value pairs:
%    None.
%
% Examples are included within the code.
%  
% See Also:
%   ieLUTDigital, ieLUTInvert
%

% History:
%    xx/xx/13       (c) Imageval Consulting, LLC 2013
%    12/06/17  jnm  Formatting

% Examples:
%{
    d = displayCreate('LCD-Apple');
    rgb = rand(10, 10, 3);
    foo = ieLUTLinear(rgb, displayGet(d, 'inverse gamma'));
    vcNewGraphWin;
    plot(foo(:), rgb(:), '.')
%}

%% Init Check inputs
%  check input parameters
if notDefined('RGB'), error('RGB value required'); end
if notDefined('invGammaTable'), invGammaTable = 2.2; end

%% Check RGB size
s = size(RGB);
nPrimaries = size(invGammaTable, 2);
if numel(RGB) == 3 || numel(RGB) == nPrimaries 
    % RGB is passed in as 3 element vector
    RGB = reshape(RGB, [1 1 length(RGB)]);
elseif ismatrix(RGB) && (size(RGB, 2) == 3 || size(RGB, 2) == nPrimaries)
    % RGB is passed in as N x 3 matrix
    RGB = reshape(RGB, [size(RGB, 1) 1 size(RGB, 2)]);
elseif ndims(RGB) > 3
    error('Unknown RGB input format');
end

%% Lookup with invert gamma table
if (numel(invGammaTable) == 1)
    % When gamma is a scalar, raise to a power
    DAC = RGB .^ invGammaTable;
else
    if size(invGammaTable, 2) == 1
        % If only one column, replicate to number of display primaries
        invGammaTable = repmat(invGammaTable, 1, size(RGB, 3));
    end
    
    % Scale the linear RGB values so that that largest value, 1 maps to the
    % row size of the gTable.
    RGB = floor(RGB * size(invGammaTable, 1)) + 1;
    RGB(RGB > size(invGammaTable, 1)) = size(invGammaTable, 1); % clip

    % Convert through the gamma table.
    DAC = zeros(size(RGB));
    for ii = 1:size(RGB, 3)
        thisTable = invGammaTable(:, ii);
        DAC(:, :, ii) = thisTable(RGB(:, :, ii));
    end 
end

%% Reshape DAC to have same size as input RGB
DAC = reshape(DAC, s);

end