function [optics, wvfP] = opticsCreate(opticsType, varargin)
% Create an optics structure
%
% Syntax:
%   [optics, wvf] = OPTICSCREATE(opticsType, [varargin])
%
% Description:
%    This function is typically called through oiCreate. The optics
%    structure is attached to the oi and manipulated by oiSet and oiGet.
%
%    The optics structure contains a variety of parameters, such as
%    f-number and focal length. There are two types of optics models:
%    diffraction limited and shift-invariant. See the discussion in
%    opticsGet for more detail.
%
%    Optics structures do not start out with a wavelength spectrum
%    structure. This information is stored in the optical image.
%
%    For diffraction-limited optics, the key parameter is the f-number.
%
%    Specifying human optics creates a shift-invariant optics structure
%    with human OTF data.
%
%    Human and general shift-invariant models can also be created by
%    specifying wavefront aberrations using Zernike polynomials. There is a
%    collection of wavefront methods to help with this (see wvfCreate,
%    wvf<TAB>). That is the method used here for 'wvf human'.
%
% Inputs:
%    opticsType - (Optional). String. Indicates the type of optics. Default
%                 is 'default'. Options listed below:
%        {'human', 'default'} - Default. Uses Marimont and Wandell
%                               (Hopkins) method. Has varargin parameter
%                               for the radius of the pupil in meters.
%        {'wvf human'}        - Uses Wavefront toolbox and Thibos data. Has
%                               varargin parameters for pupil diameter,
%                               Zernike Coefficients, wavelengths, and
%                               microns per degree.
%        {'diffraction'}      - 46 deg field of view, f number 4 optics.
%    varargin   - (Optional) Additional arguments, such as the following
%                 for a wavefront/Thibos human: (in this order)
%        pupilDiameter: Numeric. Diameter of a human pupil in millimeters.
%                       Default 3mm.
%        zCoefs:        The zernike coefficients. Default pulls from
%                       wvfLoadThibosVirtualEyes.
%        wave:          Vector. Wavelengths. Default 400:10:700.
%        umPerDegree:   Retinal parameter, microns per degree. Default 300.
%
% Outputs:
%    optics     - Struct. The created optics structure.
%    wvf        - If a wavefront optics model, the wvf struct is returned
%
% Optional key/value pairs:
%    **Needs to be filled out**
%
% Notes:
%    * TODO: Assign someone to fill in key/value section
%    * TODO: Determine if we want to implement setting transmittance freely
%      in iset
%    * TODO: Determine if planning to eventually add mouse optics.
%
% See Also:
%   oiCreate, opticsSet, opticsGet
%

% History
%	 xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    12/11/17  dhb  Allow umPerDegree to be passed for wvfhuman.
%    12/13/17  dhb  For wvfhuman, make sure that both the oi focal length
%                   and f number are set to match specified pupil diameter
%                   and umPerDegree.
%    12/20/17  dhb  Set meas pupil size for wvfHuman calcs.
%    03/13/18  jnm  Formatting
%    06/27/19  JNM  Formatting adjustments

% Examples:
%{
    % This example creates three unique optics structures based on
    % different designs/default parameters.
    optics1 = opticsCreate('diffraction');
	optics2 = opticsCreate('human');     % Marimont and Wandell
	optics3 = opticsCreate('wvf human'); % Thibos Zernike
%}

if notDefined('opticsType'), opticsType = 'default'; end

opticsType = ieParamFormat(opticsType);

switch lower(opticsType)
    case {'diffraction', 'diffractionlimited'}
        % Standard camera (cell phone) optics
        optics = opticsDiffraction;

        % For diffraction case, set lens transmittance to 1
        % Perhaps we should allow the transmittance to be set freely as in
        % ISET?  Or ...
        optics.lens = Lens;        % Human lens object
        optics.lens.density = 0;   % Pigment density set to 0

    case {'default', 'human', 'humanmw'}
        % Optics for the Marimont and Wandell human eye

        % Pupil radius in meters. Default is 3 mm
        pupilRadiusMeters = 0.0015;
        if (~isempty(varargin) && ~isempty(varargin{1}))
            pupilRadiusMeters = varargin{1};
        end

        % This creates shift-invariant optics.
        optics = opticsHuman(pupilRadiusMeters);
        optics = opticsSet(optics, 'model', 'shift invariant');
        optics = opticsSet(optics, 'name', 'human-MW');

        % Add default huamn Lens by default
        optics.lens = Lens;

    case {'wvfhuman'}
        % Optics based on mean Zernike polynomial wavefront model estimated
        % by Thibos, for 3 mm pupil. This is not actually representative of
        % any particular observer, because the mean Zernike polynomials do
        % not capture the phase information.

        % Defaults
        pupilDiameterMM = 3;
        zCoefs = wvfLoadThibosVirtualEyes(pupilDiameterMM);
        wave = 400:10:700;
        wave = wave(:);
        umPerDegree = 300;

        if (~isempty(varargin) && ~isempty(varargin{1}))
            pupilDiameterMM = varargin{1};
        end
        if (length(varargin) > 1 && ~isempty(varargin{2}))
            zCoefs = varargin{2};
        end
        if (length(varargin) > 2 && ~isempty(varargin{3}))
            wave = varargin{3};
            wave = wave(:);
        end
        if (length(varargin) > 3 && ~isempty(varargin{4}))
            umPerDegree = varargin{4};
        end

        % Create wavefront parameters. Be sure to set both measured and
        % calc pupil size.
        wvfP = wvfCreate('calc wavelengths', wave, 'zcoeffs', zCoefs, ...
            'name', sprintf('human-%d', pupilDiameterMM), ...
            'umPerDegree', umPerDegree);
        wvfP = wvfSet(wvfP, 'measured pupil size', pupilDiameterMM);
        wvfP = wvfSet(wvfP, 'calc pupil size', pupilDiameterMM);
        wvfP = wvfComputePSF(wvfP);

        optics = oiGet(wvf2oi(wvfP), 'optics');
        optics = opticsSet(optics, 'model', 'shift invariant');
        optics = opticsSet(optics, 'name', 'human-wvf');

        % Convert from pupil size and focal length to f number and focal
        % length, because that is what we can set. This implies a number of
        % mm per degree, and we back it out the other way here so that it
        % is all consistent.
        focalLengthMM = (umPerDegree * 1e-3) / (2 * tand(0.5));
        fLengthMeters = focalLengthMM * 1e-3;
        pupilRadiusMeters = (pupilDiameterMM / 2) * 1e-3;
        optics = opticsSet(optics, 'fnumber', fLengthMeters / ...
            (2 * pupilRadiusMeters));
        optics = opticsSet(optics, 'focalLength', fLengthMeters);

        % Add default Lens by default
        optics.lens = Lens;

    case 'mouse'
        % Some day might add in a default mouse optics. Here are some
        % guesses about the right parameters:
        %
        % Pupil radius in meters.
        %   Dilated pupil: 1.009mm = 0.001009m
        %   Contracted pupil: 0.178 mm
        %   (Source: From candelas to photoisomerizations in the mouse eye
        %   by rhodopsin bleaching in situ and the light-rearing dependence
        %   of the major components of the mouse ERG, Pugh, 2004)
        % We use a default value, in between: 0.59 mm.
        %         if ~isempty(varargin)
        %             pupilRadius = varargin{1};
        %             if pupilRadius > 0.001009 || pupilRadius < 0.000178
        %                 warning('Poor pupil size for the  mouse eye.')
        %             end
        %         else
        %             pupilRadius = 0.00059;  % default : 0.59 mm
        %         end
        %         % This creates a shift-invariant optics. The other
        %         % standard forms are diffraction limited.
        %         optics = opticsMouse(pupilRadius);
        %         optics = opticsSet(optics, 'model', 'shiftInvariant');
        %
        %     case {'standard(1/3-inch)', 'thirdinch'}
        %         optics = opticsThirdInch;
        %     case {'standard(1/2-inch)', 'halfinch'}
        %         optics = opticsHalfInch;
        %     case {'standard(2/3-inch)', 'twothirdinch'}
        %         optics = opticsTwoThirdInch;
        %     case {'standard(1-inch)', 'oneinch'}
        %         optics = opticsOneInch;

    otherwise
        error('Unknown optics type.');
end

% Default computational settings for the optical image
optics = opticsSet(optics, 'offAxisMethod', 'cos4th');

% Pixel vignetting is off
optics.vignetting =  0;

end

%---------------------------------------
function optics = opticsDiffraction
% Std. diffraction-limited optics w/ 46 deg. fov & fNumber of 4.
%
% Syntax:
%   optics = opticsDiffraction
%
% Description:
%    Standard diffraction limited optics with a 46-deg field of view and
%    fnumber of 4. Simple digital camera optics are like this.
%
% Inputs:
%    None.
%
% Outputs:
%    optics - Struct. The optics structure.
%
% Optional key/value pairs:
%    None.
%

optics.type = 'optics';
optics = opticsSet(optics, 'name', 'ideal (small)');
optics = opticsSet(optics, 'model', 'diffraction limited');

% Standard 1/4-inch sensor parameters
sensorDiagonal = 0.004;
FOV = 46;
fLength = inv(tan(FOV / 180 * pi) / 2 / sensorDiagonal) / 2;

optics = opticsSet(optics, 'fnumber', 4);  % focal length / diameter
optics = opticsSet(optics, 'focalLength', fLength);
optics = opticsSet(optics, 'otfMethod', 'dlmtf');

end

%---------------------------------------
function optics = opticsHuman(pupilRadius)
% Using shift-invariant, place estimated human OTF with humanOTF
%
% Syntax:
%   optics = opticsHuman([pupilRadius])
%
% Description:
%    Use the shift-invariant method place the estimated human OTF, using
%    humanOTF from Marimont and Wandell, in the OTF fields.
%
%    The frequency support is calculated in cyc/deg but stored in units of
%    cyc/mm.
%
%    In the wvf code, use we 300 microns/deg as the conversion factor. This
%    value corresponds to a focal length of 17.188 mm or the human eye.
%    Here we use 17 mm, and a 3 mm pupil.
%
% Inputs:
%    pupilRadius - (Optional) Numeric. Radius of pupil in meters. Default
%                  is 0.0015m.
%
% Outputs:
%    optics      - Struct. The created optics structure.
%
% Optional key/value pairs:
%    None.
%

% The pupil radius is specified in meters.
if notDefined('pupilRadius'), pupilRadius = 0.0015; end

% Human focal length is ~17 mm. This corresponds to 296.71 um per degree.
% Elsewhere we use 300 um per degree, but 17mm is what we've had here and
% is what the Marimont and Wandell optics is based on so we keep that here.
focalLengthMM = 17.1883;
focalLengthMM = 17;
fLengthMeters = focalLengthMM * 1e-3;

% Calculate umPerDegree. We don't use this here, but this was on the way to
% understanding what various numbers we have hard coded into various places
% of the code.
mmPerDegree = 2 * focalLengthMM * tand(0.5);
umPerDegree = mmPerDegree * 1e3;

% Start setting up the optics structure.
optics.type = 'optics';
optics.name = 'human';
optics = opticsSet(optics, 'model', 'shiftInvariant');

% Convert from pupil size and focal length to f number and focal length,
% because that is what we can set. This implies a number of mm per degree.
optics = opticsSet(optics, 'fnumber', fLengthMeters / (2 * pupilRadius));
optics = opticsSet(optics, 'focalLength', fLengthMeters);

% Specify method for optics
optics = opticsSet(optics, 'otfMethod', 'humanOTF');

% Compute the OTF and store it. We use a default pupil radius, dioptric
% power, and so forth.
dioptricPower = 1 / fLengthMeters;  % About 60 diopters

% We used to assign the same wave as in the current scene to optics, if the
% wave was not yet assigned.
wave = opticsGet(optics, 'wave');

% Decide whether to use the legacy frequency support (Nyquist:  60 c/deg)
% or 120 c/deg, which results in a more focused PSF
legacyFrequencySupport = true;
if legacyFrequencySupport
    % Fsupport used to be [], which defaults to 60 c/deg
    fSupport = [];
else
    % Up to 120 c/deg, with 240 samples to get better PSF
    maxF = 120;
    N = 240;
    fList = unitFrequencyList(N);
    fList = fList * maxF;
    [X, Y] = meshgrid(fList, fList);
    fSupport(:, :, 1) = X;
    fSupport(:, :, 2) = Y;
end

% The human optics are an SI case, and we store the OTF at this point.
[OTF2D, fSupport] = humanOTF(pupilRadius, dioptricPower, fSupport, wave);

optics = opticsSet(optics, 'otfData', OTF2D);
umPerDegreeForSupport = umPerDegree;

fSupport = fSupport * (1 / (umPerDegreeForSupport * 1e-3));
fx = fSupport(1, :, 1);
fy = fSupport(:, 1, 2);
optics = opticsSet(optics, 'otffx', fx(:)');
optics = opticsSet(optics, 'otffy', fy(:)');
optics = opticsSet(optics, 'otfWave', wave);

end