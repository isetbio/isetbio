function oi = oiCreate(oiType, varargin)
% Create an optical image structure
%
% Syntax:
%   oi = oiCreate(oiType, [varargin])
%
% Description:
%    The optical image represents the spectral irradiance at the sensor.
%    The irradiance is computed from the scene radiance, using the
%    information in the optics structure that is attached to this oi.
%
%    The oi spectrum is normally inherited from the scene. To specify a
%    spectrum for the optical image use
%        oi = oiCreate('default');
%        oi = initDefaultSpectrum('hyperspectral');
%
%    Lens transmittance
%       The human models include a lens object within the optics structure.
%       The lens object represents information about the lens pigment.
%
%       The diffraction model does not include the lens object; instead, it
%       includes a transmittance structure with a 'wave' and 'scale' slot.
%       These slots characterize the lens spectral transmittance. The
%       values in the transmittance structure are interpolated as needed
%       during the calculation.
%
%    There are examples contained in the code below. To access, type 'edit
%    oiCreate.m' into the Command Window.
%
% Inputs:
%    oiType   - (Optional) String. The type of oi structure. The two
%               options marked with * are typically used for testing and
%               not called for experiments. Options are:
%       {'human'}      - (Default) Human shift-invariant optics based on
%                        Marimont and Wandell (1994, JOSA)
%       {'wvf human'}  - Human shift-invariant optics based on mean
%                        wavefront abberration from Thibos et al. (2009,
%                        Ophthalmic & Physiological Optics). Optional
%                        parameters can be passed for this case (see below)
%      {'diffraction'} - Diffraction limited optics, f / 4, contains no
%                        diffuser or data.
%      {'uniformd65'}* - Turns off offaxis to create an uniform D65 image
%      {'uniformee'}*  - Turns off offaxis to create an uniform equal
%                        energy image.
%    varargin - (Optional) Additional parameters based on the oiType. They
%               are listed with their type in the required order and
%               default values below:
%      'wvf human' - A human optical image structurem such as contained in
%                    oiCreate('wvf human', pupilMM, zCoefs, wave, umPerDeg)
%                     pupilMM:     Numerical. Pupil Size. Default 3mm.
%                     zCoefs:      () z Coefficients Default calls the
%                                  function wvfLoadThibosVirtualEyes.
%                     wave:        Vector. Wavelengths. Default 400:10:700.
%                     umPerDegree: Scalar. Microns per Degree. Default 300.
%
% Outputs:
%    oi       - Struct. The created optical image structure
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * TODO: Determine if still need to set bit depth to single precision.
%
% See Also:
%   sceneCreate, opticsCreate
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    03/06/18  jnm  Formatting
%    06/25/19  JNM  Minor formatting adjustments

% Examples:
%{
    % Marimont/Wandell model
    oi = oiCreate('human');
%}
%{
    oi = oiCreate('diffraction');
    oiPlot(oi, 'lsf wavelength');
%}
%{
    % Mean Thibos et al. data. Not too blurry.
    oi = oiCreate('wvf human');
    pupilMM = 2;
    oi = oiCreate('wvf human', pupilMM)
    oiPlot(oi, 'lsf wavelength');
    title(sprintf('%d mm', pupilMM))
%}
%{
    % Bigger pupil, much blurrier. You can barely see outside of 520-580 nm
    oi = oiCreate('wvf human');
    pupilMM = 5;
    oi = oiCreate('wvf human', pupilMM)
    oiPlot(oi, 'lsf wavelength');
    title(sprintf('%d mm', pupilMM))
%}

%%
if notDefined('oiType'), oiType = 'human'; end

% Default is to use the diffraction limited calculation
oi.type = 'opticalimage';
oi.name = vcNewObjectName('opticalimage');  % Get a fresh name
% Set depth to single precision. Perhaps no longer needed?
oi = oiSet(oi, 'bit depth', 32);

oiType = ieParamFormat(oiType);
switch oiType
    case {'default', 'human', 'mwhuman'}
        % Marimont and Wandell optics, which is a simple shift-invariant
        % but wavelength-dependent model. This is a little faster than the
        % wvf human, so we made it the default. They differ a little.
        %
        % oi = oiCreate('human');
        oi = oiCreate('diffraction limited');
        oi = oiSet(oi, 'diffuser method', 'skip');
        oi = oiSet(oi, 'consistency', 1);
        oi = oiSet(oi, 'optics', opticsCreate('human'));
        oi = oiSet(oi, 'name', 'human-MW');
        oi = oiSet(oi, 'lens', Lens('wave', oiGet(oi, 'optics wave')));

        if checkfields(oi.optics, 'transmittance')
            oi.optics = rmfield(oi.optics, 'transmittance');
        end

    case {'wvfhuman', 'shiftinvariant'}
        % A human lens specified using the WVF toolbox method
        % oi = oiCreate('wvf human', pupilMM, zCoefs, wave)

        oi = oiCreate('diffraction limited');
        oi = oiSet(oi, 'diffuser method', 'skip');
        oi = oiSet(oi, 'consistency', 1);

        oi = oiSet(oi, 'optics', opticsCreate('wvf human', varargin{:}));
        oi = oiSet(oi, 'name', 'human-WVF');
        oi = oiSet(oi, 'lens', Lens('wave', oiGet(oi, 'optics wave')));

        if checkfields(oi.optics, 'transmittance')
            oi.optics = rmfield(oi.optics, 'transmittance');
        end

    case {'diffractionlimited', 'diffraction'}
        % Default optics is f# = 4, diffraction limited
        optics = opticsCreate('diffraction limited');
        oi = oiSet(oi, 'optics', optics);
        oi = oiSet(oi, 'name', 'diffraction');

        % Set up the default glass diffuser with a 2 micron blur circle,
        % but skipped
        oi = oiSet(oi, 'diffuser method', 'skip');
        oi = oiSet(oi, 'diffuser blur', 2 * 10 ^ -6);
        oi = oiSet(oi, 'consistency', 1);

        if checkfields(oi.optics, 'lens')
            oi.optics = rmfield(oi.optics, 'lens');
            oi.optics.transmittance.wave = (370:730)';
            oi.optics.transmittance.scale = ones(length(370:730), 1);
        end

    % The following two are used for testing. They are not really
    % models.
    case {'uniformd65'}
        % No lens, just transmittance structure (diffraction type) Uniform,
        % D65 optical image. No cos4th falloff, huge field of view (120
        % deg). Used in lux-sec SNR testing and scripting
        oi = oiCreateUniformD65;

    case {'uniformee'}
        % No lens, just transmittance structure (diffraction type)
        % Uniform, equal energy optical image. No cos4th falloff, huge
        % field of view (120 deg).
        oi = oiCreateUniformEE;

    otherwise
        error('Unknown oiType');
end

end

function oi = oiCreateUniformD65
% Create a uniform D65 image w/ large FoV.
%
% Syntax:
%   oi = oiCreateUniformD65;
%
% Description:
%    Create a uniform, D65 image with a very large field of view. The
%    optical image is created without any cos4th fall off so it can be used
%    for lux-sec SNR testing.
%
% Inputs:
%    None.
%
% Outputs:
%    oi - The created optical image structure
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * [Note: XXX - This does not yet extend in the IR, but it should. See
%      notes in sceneCreate.]
%

scene = sceneCreate('uniform d65');
scene = sceneSet(scene, 'hfov', 120);

oi = oiCreate('diffraction');
oi = oiSet(oi, 'optics fnumber', 1e-3);   % Basically perfect optics
oi = oiCompute(scene, oi);
% vcAddObject(oi);
% oiWindow;

end

function oi = oiCreateUniformEE
% Create a uniform equal energy image with a large FoV
%
% Syntax:
%   oi = oiCreateUniformEE;
%
% Description:
%    Create a uniform, equal energy image with a very large field of view.
%    The optical image is created without any cos4th fall off so it can be
%    used for lux-sec SNR testing.
%
% Inputs:
%    None.
%
% Outputs:
%    oi - The created optical image structure
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * [Note: XXX - This does not yet extend in the IR, but it should. See
%      notes in sceneCreate.]
%

scene = sceneCreate('uniform ee');
scene = sceneSet(scene, 'hfov', 120);

oi = oiCreate('diffraction');
oi = oiSet(oi, 'optics fnumber', 1e-3);   % Basically perfect optics
oi = oiCompute(scene, oi);
% vcAddObject(oi);
% oiWindow;

end
