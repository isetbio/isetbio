function LMS = lms2lmsDichromat(LMS, cbType, method, varargin)
%% Interpolate missing cone values for colorblind in cone color space (LMS)
%
% Syntax:
%   LMS = lms2lmsDichromat(LMS, cbType, method, varargin)
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
%       proposed by Jiang, Joyce and Wandell, 2015
%     * returning a zero for the missing cone type  
%
% Inputs:
%    LMS    - an image in cone space (Stockman LMS Space)
%    cbType - Type of colorblindness, can be choosen from
%      {0, 'Trichromat'}             - Trichromatic observer
%      {1, 'Protan', 'Protanopia'}   - Protanopia observer, missing L cones
%      {2, 'Deutan', 'Deuteranopia'} - Deuteranope, missing M cones
%      {3, 'Tritan', 'Tritanopia'}   - Tritanope, missing S cones
%    method - Algorithm to be used to interpolate the missing cone values
%      'Brettel'  - Using Bettel's algorithm (Default)
%      'Linear'   - Using Linear Interpolation
%      'Constant' - leave missing cone values as a constant (default 0)
%
% Optional Key/Value Pairs:
%    varargin - More parameters input for different methods
%      varargin{1} - white LMS values for Brettel mehtod
%                    extrapolate values for Constant method
%
% Outputs:
%    LMS - LMS values in Stockman LMS space
%
% See Also:
%    xyz2lms
%
% (HJ) ISETBIO TEAM, 2015

% Examples:
%{
   LMS = lms2lmsDichromat(LMS, colorBlindType, 'linear');
%}
%% Init and Check Inputs
if notDefined('LMS'), error('LMS image required'); end
if notDefined('cbType'), cbType = 0; end
if notDefined('method'), method = 'Brettel'; end

if ischar(cbType), cbType = ieParamFormat(cbType); end
method = ieParamFormat(method);


%% Interpolate missing cone
switch method
    case 'brettel'
        if isempty(varargin), error('white LMS values required'); end
        LMS = cbBrettel(LMS, cbType, varargin{1});
    case 'linear'
        LMS = cbLinear(LMS, cbType);
    case 'constant'
        if isempty(varargin), val = 0; else val = varargin{1}; end
        LMS = cbConstant(LMS, cbType, val);
end
end

%% Interpolation with Brettel's Method
function LMS = cbBrettel(LMS, cbType, whiteLMS)
% Interpolation with Brettel's Method
%
% Syntax:
%   LMS = cbBrettel(LMS, cbType, whiteLMS)
%
% Description:
%    Perform your interpolation using Brettel's Method
%
% Inputs:
%    LMS      - LMS Color space values
%    cbType   - Type of Color-blindness
%    whiteLMS - The LMS White Point
%
% Outputs:
%    LMS      - LMS Color space values after manipulation to account for
%               the specified type of color blindness
%
% Notes:
%    * [Note: JNM - Unsure if the definitions I have supplied for the
%       inputs and output are sufficient]

% Examples:
%{
   LMS = cbBrettel(LMS, cbType, whiteLMS)
%}

% check inputs
if notDefined('LMS'), error('LMS image required'); end
if notDefined('cbType'), cbType = 0; end
if notDefined('whiteLMS'), error('white LMS values required'); end
if ischar(cbType), cbType = ieParamFormat(cbType); end

% define the anchor points
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
    case {0, 'trichromat'} % Trichromatic Observer, do nothing
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
        
        % Interpolate missing L values for protonate
        L = LMS(:,:,1); M = LMS(:,:,2); S = LMS(:,:,3);
        lst = ((S ./ M) < inflection);
        L(lst)  = -(b1*M(lst)  + c1*S(lst))  / a1;
        L(~lst) = -(b2*M(~lst) + c2*S(~lst)) / a2;
        LMS(:,:,1) = L;
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
        L = LMS(:,:,1); M = LMS(:,:,2); S = LMS(:,:,3);
        lst = ((S ./ L) < inflection);
        M(lst)  = -(a1*L(lst)  + c1*S(lst)) / b1;
        M(~lst) = -(a2*L(~lst) + c2*S(~lst))/ b2;
        LMS(:,:,2) = M;
    case {3, 'tritan', 'tritanope', 'tritanopia'}
        a1 = whiteLMS(2) * anchor(12) - whiteLMS(3) * anchor(11);
        b1 = whiteLMS(3) * anchor(10)  - whiteLMS(1) * anchor(12);
        c1 = whiteLMS(1) * anchor(11) - whiteLMS(2) * anchor(10);
        
        % Greater than the inflection
        a2 = whiteLMS(2) * anchor(6)  - whiteLMS(3) * anchor(5);
        b2 = whiteLMS(3) * anchor(4)  - whiteLMS(1) * anchor(6);
        c2 = whiteLMS(1) * anchor(5)  - whiteLMS(2) * anchor(4);
        
        % Inflection point
        inflection = (whiteLMS(2) / whiteLMS(1));
        
        % Interpolate missing M values for tritanope
        L = LMS(:,:,1); M = LMS(:,:,2); S = LMS(:,:,3);
        lst = ((M ./ L) < inflection);
        S(lst)  = -(a1*L(lst)  + b1*M(lst)) / c1;
        S(~lst) = -(a2*L(~lst) + b2*M(~lst))/ c2;
        LMS(:,:,3) = S;
    otherwise
        error('unknown color blind type');
end
end

%% Interpolation with Linear Method
function LMS = cbLinear(LMS, cbType)
% Interpolation with Linear Method
%
% Syntax:
%   LMS = cbLinear(LMS, cbType)
%
% Description:
%    Interpolation with Linear Method
%
% Inputs:
%    LMS      - LMS Color space values
%    cbType   - Type of Color-blindness
%
% Outputs:
%    LMS      - LMS Color space values after manipulation to account for
%               the specified type of color blindness
%
% Notes:
%    * [Note: JNM - Unsure if the definitions I have supplied for the
%       inputs and output are sufficient]

% Examples:
%{
   LMS = cbLinear(LMS, cbType)
%}

% check inputs
if notDefined('LMS'), error('LMS image required'); end
if notDefined('cbType'), cbType = 0; end
if ischar(cbType), cbType = ieParamFormat(cbType); end

% Interpolate missing cone values
switch cbType
    case {0, 'trichromat'} % Trichromatic Observer, do nothing
    case {1, 'protan', 'protanopia', 'protanope'}
        LMS(:,:,1) = 1.3855*LMS(:,:,2) - 0.285*LMS(:,:,3);
    case {2, 'deutan', 'deuteran', 'deuteranope', 'deuteranopia'}
        LMS(:,:,2) = 0.6949*LMS(:,:,1) + 0.2614*LMS(:,:,3);
    case {3, 'tritan', 'tritanope', 'tritanopia'}
        LMS(:,:,3) = -0.9623*LMS(:,:,1) + 1.7595*LMS(:,:,2);
    otherwise
        error('unknown color blind type');
end
end

%% Interpolation with Constants
function LMS = cbConstant(LMS, cbType, val)
% Interpolation with Constants
%
% Syntax:
%   LMS = cbConstant(LMS, cbType, val)
%
% Description:
%    Interpolation with a Constant
%
% Inputs:
%    LMS      - LMS Color space values
%    cbType   - Type of Color-blindness
%    val      - The requisite constant
%
% Outputs:
%    LMS      - LMS Color space values after manipulation to account for
%               the specified type of color blindness
%
% Notes:
%    * [Note: JNM - Unsure if the definitions I have supplied for the
%       inputs and output are sufficient]

% Examples:
%{
   LMS = cbConstant(LMS, cbType, val)
%}

% check inputs
if notDefined('LMS'), error('LMS image required'); end
if notDefined('cbType'), cbType = 0; end
if notDefined('val'), val = 0; end
if ischar(cbType), cbType = ieParamFormat(cbType); end

% Interpolate missing cone values
switch cbType
    case {0, 'trichromat'} % Trichromatic Observer, do nothing
    case {1, 'protan', 'protanopia', 'protanope'}
        LMS(:,:,1) = val;
    case {2, 'deutan', 'deuteran', 'deuteranope', 'deuteranopia'}
        LMS(:,:,2) = val;
    case {3, 'tritan', 'tritanope', 'tritanopia'}
        LMS(:,:,3) = val;
    otherwise
        error('unknown color blind type');
end
end
