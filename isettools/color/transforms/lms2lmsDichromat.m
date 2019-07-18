function LMS = lms2lmsDichromat(LMS, cbType, method, varargin)
% Interpolate missing cone values for colorblind in cone color space (LMS)
%
% Syntax:
%   LMS = lms2lmsDichromat(LMS, [cbType], [method], [varargin])
%
% Description:
%    Interpolate missing cone values for colorblind in cone color space
%
%    The estimate for the missing cone type for colorblind is done by
%     * interpolating the missing cone using the algorithm in Brettel,
%       Vienot and Mollon JOSA 14/10 1997. The idea in that paper is that
%       the preserved cones are preserved. The missing cone is assigned a
%       value that is a piece-wise linear transform of the preserved cones
%     * interpolation the missing cone using the linear interpolation
%       proposed by Jiang, Farrell and Wandell, 2015
%     * returning a constant for the missing cone type
%
%    This function contains examples of usage inline. To access these, type
%    'edit lms2lmsDichromat.m' into the Command Window.
%
% Inputs:
%    LMS      - Matrix. An image in cone space (Stockman LMS Space).
%    cbType   - (Optional) Scalar/String. A type of colorblindness, can be
%               choosen from the following:
%      {0, 'trichromat'}             - Default. Trichromatic observer.
%      {1, 'protan', 'protanopia'}   - Protanopia observer, missing L cones
%      {2, 'deutan', 'deuteranopia'} - Deuteranope, missing M cones
%      {3, 'tritan', 'tritanopia'}   - Tritanope, missing S cones
%    method   - (Optional) String. Which algorithm to be used to
%               interpolate the missing cone values. Options (and their
%               varargin information) include:
%        brettel: Default. Use Bettel's algorithm. Provide numeric value to
%                 varargin{1} of a white XYZ value.
%        linear: Use linear interpolation algorithm
%        constant: Leave missing cone values as a constant. Default 0.
%    varargin - (Optional) VARIES. Extra arguments for specific methods.
%               See the method entries above for more information.
%
% Outputs:
%    LMS    - Matrix. The LMS values in Stockman LMS space.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * [NOTE - DHB: I made the example run, but did not actually check that
%      the output is correct.]
%    * [NOTE - DHB: Not clear if it is 2 degree or 10 degree XYZ/Stockman.
%      There are various constants hard coded into the transformations, so
%      to determine which it is one would need to recreate where those
%      constants came from.]
%
% See Also:
%   xyz2lms
%

% History:
%    xx/xx/15  HJ   ISETBIO TEAM
%    11/01/17  jnm  Comments, formatting & fix example
%    11/17/17  jnm  Formatting & fix note
%    07/16/19  JNM  Formatting update

% Examples:
%{
    scene = sceneCreate;
    imgLMS = sceneGet(scene, 'lms');
    LMS = lms2lmsDichromat(imgLMS, 0, 'linear');
    xyz2lmsM = colorTransformMatrix('xyz2lms');
    whiteXYZ = sceneGet(scene,'illuminant xyz');
    whiteLMS = whiteXYZ(:)' * xyz2lmsM;
    LMS = lms2lmsDichromat(imgLMS, 'deutan', 'brettel', whiteLMS);
%}

%% Init and check inputs
if notDefined('LMS'), error('LMS image required'); end
if notDefined('cbType'), cbType = 0; end
if notDefined('method'), method = 'Brettel'; end

if ischar(cbType), cbType = ieParamFormat(cbType); end
method = ieParamFormat(method);

%% Interpolate missing cone
switch method
    case 'brettel'
        if isempty(varargin), error('White point LMS values required'); end
        LMS = cbBrettel(LMS, cbType, varargin{1});
    case 'linear'
        LMS = cbLinear(LMS, cbType);
    case 'constant'
        if isempty(varargin)
            val = 0;
        else
            val = varargin{1};
        end
        LMS = cbConstant(LMS, cbType, val);
end
end

%% Interpolation with Brettel's method
function LMS = cbBrettel(LMS, cbType, whiteLMS)
% Interpolation with Brettel's method
%
% Syntax:
%   LMS = cbBrettel(LMS, [cbType], whiteLMS)
%
% Description:
%    Perform the interpolation using Brettel's Method
%
% Inputs:
%    LMS      - Matrix. LMS color space values
%    cbType   - (Optional) Scalar/String. Type of color-blindness. See the
%               variable of the same name above in lms2lmsDichromat.
%               Default 0.
%    whiteLMS - Vector. The LMS white point.
%
% Outputs:
%    LMS      - Matrix. LMS Color space values after manipulation to
%               account for the specified type of color blindness
%
% Optional key/value pairs:
%    None.
%

% Check inputs
if notDefined('LMS'), error('LMS image required'); end
if notDefined('cbType'), cbType = 0; end
if notDefined('whiteLMS'), error('white LMS values required'); end
if ischar(cbType), cbType = ieParamFormat(cbType); end

% Define the anchor points
% For protan, the anchor wavelength are 475 nm and 575 nm. For deuteranope,
% the anchor wavelength are 475 nm and 575 nm. For tritanope, the anchor
% wavelength are 485 nm and 660 nm.
%
% The anchor values are roughly
%   anchor = ieReadSpectra('stockman', [475 485 575 660]);
anchor = [0.1193 0.2123 0.5309
          0.1620 0.2709 0.2894
          0.9881 0.7266 0.0004
          0.0941 0.0077 0.0]';
% anchor = [0.1420    0.2413    0.4964
%           0.1888    0.3016    0.2650
%           0.9716    0.6824    0.0003
%           0.0807    0.0063    0.0000]';

% Interpolate missing cone values
switch cbType
    case {0, 'trichromat'} % Trichromatic observer, do nothing
    case {1, 'protan', 'protanopia', 'protanope'}
        a1 = whiteLMS(2) * anchor(9) - whiteLMS(3) * anchor(8);
        b1 = whiteLMS(3) * anchor(7) - whiteLMS(1) * anchor(9);
        c1 = whiteLMS(1) * anchor(8) - whiteLMS(2) * anchor(7);
        % Greater than inflection
        a2 = whiteLMS(2) * anchor(3) - whiteLMS(3) * anchor(2);
        b2 = whiteLMS(3) * anchor(1) - whiteLMS(1) * anchor(3);
        c2 = whiteLMS(1) * anchor(2) - whiteLMS(2) * anchor(1);

        % Divides the space
        inflection = whiteLMS(3) / whiteLMS(2);

        % Interpolate missing L values for protonape
        L = LMS(:, :, 1);
        M = LMS(:, :, 2);
        S = LMS(:, :, 3);
        lst = ((S ./ M) < inflection);
        L(lst)  = -(b1 * M(lst)  + c1 * S(lst))  / a1;
        L(~lst) = -(b2 * M(~lst) + c2 * S(~lst)) / a2;
        LMS(:, :, 1) = L;
    case {2, 'deutan', 'deuteran', 'deuteranope', 'deuteranopia'}
        a1 = whiteLMS(2) * anchor(9) - whiteLMS(3) * anchor(8);
        b1 = whiteLMS(3) * anchor(7) - whiteLMS(1) * anchor(9);
        c1 = whiteLMS(1) * anchor(8) - whiteLMS(2) * anchor(7);
        % Greater than inflection
        a2 = whiteLMS(2) * anchor(3) - whiteLMS(3) * anchor(2);
        b2 = whiteLMS(3) * anchor(1) - whiteLMS(1) * anchor(3);
        c2 = whiteLMS(1) * anchor(2) - whiteLMS(2) * anchor(1);

        inflection = whiteLMS(3) / whiteLMS(1);

        % Interpolate missing M values for deuteranope
        L = LMS(:, :, 1);
        M = LMS(:, :, 2);
        S = LMS(:, :, 3);
        lst = ((S ./ L) < inflection);
        M(lst)  = -(a1 * L(lst)  + c1 * S(lst)) / b1;
        M(~lst) = -(a2 * L(~lst) + c2 * S(~lst)) / b2;
        LMS(:, :, 2) = M;
    case {3, 'tritan', 'tritanope', 'tritanopia'}
        a1 = whiteLMS(2) * anchor(12) - whiteLMS(3) * anchor(11);
        b1 = whiteLMS(3) * anchor(10) - whiteLMS(1) * anchor(12);
        c1 = whiteLMS(1) * anchor(11) - whiteLMS(2) * anchor(10);

        % Greater than the inflection
        a2 = whiteLMS(2) * anchor(6)  - whiteLMS(3) * anchor(5);
        b2 = whiteLMS(3) * anchor(4)  - whiteLMS(1) * anchor(6);
        c2 = whiteLMS(1) * anchor(5)  - whiteLMS(2) * anchor(4);

        % Inflection point
        inflection = (whiteLMS(2) / whiteLMS(1));

        % Interpolate missing M values for tritanope
        L = LMS(:, :, 1);
        M = LMS(:, :, 2);
        S = LMS(:, :, 3);
        lst = (M ./ L) < inflection;
        S(lst)  = -(a1 * L(lst)  + b1 * M(lst)) / c1;
        S(~lst) = -(a2 * L(~lst) + b2 * M(~lst)) / c2;
        LMS(:, :, 3) = S;
    otherwise
        error('unknown color blind type');
end
end

%% Interpolation with linear method
function LMS = cbLinear(LMS, cbType)
% Interpolation with linear method
%
% Syntax:
%   LMS = cbLinear(LMS, [cbType])
%
% Description:
%    Interpolation with linear method
%
% Inputs:
%    LMS    - Matrix. The LMS color space values.
%    cbType - (Optional) Scalar/Numeric. Type of color-blindness. Default
%             0. See the variable of the same name above in
%             lms2lmsDichromat for more information.
%
% Outputs:
%    LMS    - Matrix. The LMS color space values after manipulation to
%             account for the specified type of color blindness.
%
% Optional key/value pairs:
%    None.
%

% Check inputs
if notDefined('LMS'), error('LMS image required'); end
if notDefined('cbType'), cbType = 0; end
if ischar(cbType), cbType = ieParamFormat(cbType); end

% Interpolate missing cone values
switch cbType
    case {0, 'trichromat'} % Trichromatic Observer, do nothing
    case {1, 'protan', 'protanopia', 'protanope'}
        LMS(:, :, 1) = 1.3855 * LMS(:, :, 2) - 0.285 * LMS(:, :, 3);
    case {2, 'deutan', 'deuteran', 'deuteranope', 'deuteranopia'}
        LMS(:, :, 2) = 0.6949 * LMS(:, :, 1) + 0.2614 * LMS(:, :, 3);
    case {3, 'tritan', 'tritanope', 'tritanopia'}
        LMS(:, :, 3) = -0.9623 * LMS(:, :, 1) + 1.7595 * LMS(:, :, 2);
    otherwise
        error('unknown color blind type');
end
end

%% Interpolation with constant
function LMS = cbConstant(LMS, cbType, val)
% Interpolation with constant
%
% Syntax:
%   LMS = cbConstant(LMS, [cbType], val)
%
% Description:
%    Interpolation with a constant value. Just use that value for the
%    missing cone type.
%
% Inputs:
%    LMS    - Matrix. The LMS color space values.
%    cbType - (Optional) Scalar/String. The type of color-blindness.
%             Default 0. See the variable of the same name above in
%             lms2lmsDichromat for more information.
%    val    - Numeric. A scalar numeric of the requisite constant.
%
% Outputs:
%    LMS    - Matrix. The LMS color space values after manipulation to
%             account for the specified type of color blindness
%
% Optional key/value pairs:
%    None.
%

% check inputs
if notDefined('LMS'), error('LMS image required'); end
if notDefined('cbType'), cbType = 0; end
if notDefined('val'), val = 0; end
if ischar(cbType), cbType = ieParamFormat(cbType); end

% Interpolate missing cone values
switch cbType
    case {0, 'trichromat'} % Trichromatic Observer, do nothing
    case {1, 'protan', 'protanopia', 'protanope'}
        LMS(:, :, 1) = val;
    case {2, 'deutan', 'deuteran', 'deuteranope', 'deuteranopia'}
        LMS(:, :, 2) = val;
    case {3, 'tritan', 'tritanope', 'tritanopia'}
        LMS(:, :, 3) = val;
    otherwise
        error('unknown color blind type');
end
end
