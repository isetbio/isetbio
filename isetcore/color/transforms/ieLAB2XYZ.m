function xyz = ieLAB2XYZ(lab, whitepoint, useOldCode, labexp)
% Convert CIE LAB values to CIE XYZ values
%
% Syntax:
%   xyz = ieLAB2XYZ(lab, whitepoint, [useOldCode], [labexp])
%
% Description:
%    Converts CIE L*a*b* coordinates to CIE XYZ coordinates. We will use
%    the makecform routine from the Matlab image processing toolbox for the
%    converison; if the toolbox/routine is not available, we will revert to
%    the older version of the code.
%
%    Either the 2 or 10 degree XYZ values may be used, depending on field
%    size. If the size is less than 4 degrees, then the 2-deg fundamentals
%    are recommended. If more than 4 degrees, then the 10-deg
%    fundamentals.
%
%    This function contains examples of usage inline. To access these, type
%    'edit ieLAB2XYZ.m' into the Command Window.
%
% Inputs:
%    lab        - Matrix. The LAB image; can either be in XW or RGB format.
%    whitepoint - Vector. A 3-vector of the xyz values of the white point.
%    useOldCode - (Optional) Boolean. A numeric boolean to indicate whether
%                 or not to use old code. The options are either 0 to use
%                 Matalb's routines, (1) otherwise. Default 0.
%    labexp     - (Optional) Numeric. used by old code; the exponent used
%                 in the CIELAB formula Default of 3 represents the cube
%                 root as used in standard CIELAB. If specified, use the
%                 number as exponent.
%
% Outputs:
%    xyz        - Matrix. The CIE XYZ values.
%
% Optional key/value pairs:
%    None.
%
% See Also:
%   ieXYZ2LAB
%

% History
%    08/18/15  dhb  Change conditional on exist of makecform, works for
%                   p-code too.
%              dhb  Always define labexp, since makecform may not exist.
%              dhb  Change "exp"->"labexp" to avoid clobbering function exp
%    10/25/17  jnm  Comments & Formatting
%    11/01/17  jnm  Fixed final reference to exp instead of labexp
%    11/17/17  jnm  Formatting & fix example (useOldCode was not
%                   instantiated, replaced with a 0)
%                   Copyright ImagEval Consultants, LLC, 2009.
%    01/20/18  dhb  Fix example so it works.
%    07/17/19  JNM  Formatting update

% Examples:
%{
    dataXYZ = [100 100 100 ; 50 100 75; 80 90 110];
    whiteXYZ = [100 100 100];
    lab = ieXYZ2LAB(dataXYZ, whiteXYZ)
    xyz = ieLAB2XYZ(lab, whiteXYZ)
    xyz = ieLAB2XYZ(lab, whiteXYZ, true)
%}

if notDefined('lab'), error('No data.'); end
if notDefined('whitepoint')
    error('A whitepoint is required for conversion to CIELAB (1976).');
end
if notDefined('useOldCode'), useOldCode = 0; end
if notDefined('labexp'), labexp = 3; end

if exist('makecform', 'file') && ~useOldCode
    % Which version of LAB is this for? 1976.
    % We are worried about the white point.
    cform = makecform('lab2xyz', 'WhitePoint', whitepoint(:)');
    xyz = applycform(lab, cform);
    return;
else
    if length(whitepoint) ~= 3
        error('White point is not a three-vector');
    else
        Xn = whitepoint(1);
        Yn = whitepoint(2);
        Zn = whitepoint(3);
    end

    % We will always work in XW format. If input is in RGB format, we
    % reshape it
    if ndims(lab) == 3
        [r, c, ~] = size(lab);
        lab = RGB2XWFormat(lab);
    end

    % Usual formula for Lstar. (y = Y/Yn)
    fy = (lab(:, 1) + 16) / 116;
    y = fy .^ labexp;

    % Find out cases where (Y/Yn) is too small and use other formula
    % Y / Yn = 0.008856 correspond to L = 7.9996
    yy = find(lab(:, 1) <= 7.9996);
    y(yy) = lab(yy, 1) / 903.3;
    fy(yy) = 7.787 * y(yy) + 16 / 116;

    % find out fx, fz
    fx = lab(:, 2) / 500 + fy;
    fz = fy - lab(:, 3) / 200;

    % find out x = X / Xn, z = Z / Zn
    % when (X / Xn) < 0.008856, fx < 0.206893
    % when (Z / Zn) < 0.008856, fz < 0.206893
    xx = find(fx < .206893);
    zz = find(fz < .206893);
    x = fx .^ labexp;
    z = fz .^ labexp;
    x(xx) = (fx(xx) - 16 / 116) / 7.787;
    z(zz) = (fz(zz) - 16 / 116) / 7.787;

    xyz = [x * Xn, y * Yn, z * Zn];

    % Return XYZ in appropriate shape
    if ndims(xyz) == 3, xyz = XW2RGBFormat(xyz, r, c); end
end
end
