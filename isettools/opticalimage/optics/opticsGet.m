function val = opticsGet(optics, parm, varargin)
% Get optics parameters (some values are stored but many are computed)
%
% Syntax:
%   val = opticsGet(optics, parm, [varargin])
%
% Description:
%    The optics parameters are stored in two different groups of fields.
%
%    There are two optics models implemented. The method can be selected by
%    the popup menu in the optics window or programatically via
%
%       oiSet(oi, 'optics model', <model type>);
%
%    (1) By default, we  use a general shift-invariant calculation
%    based on a numerically-defined OTF that is wavelength-dependent
%    (but shift-invariant), stored in the optics.OTF structure. When
%    using this method, the user can supply the optics structure
%    containing the OTF and other parameters (focal length, aperture,
%    and so forth). Examples are stored in data/optics directory. It
%    is also possible to create shift invariant OTF data from
%    wavefront aberrations specified as Zernike polynomials. In this
%    case, the optical image is computed using the function
%    opticsSICompute. To set this method use oi = oiSet(oi, 'optics
%    model', 'shift invariant').
%
%    (2) We use include a diffraction limited optics with the f-number
%    and focal length determined from the window interface. Other
%    parameters are derived from these two values. In this case, the
%    optical image is computed using oiCompute. To set this method,
%    use opticsSet(optics, 'model', 'diffractionLimited'); In this
%    case, the optical image is computed using opticsDLCompute. To set
%    this method use oi = oiSet(oi, 'optics model', 'diffraction
%    limited').
%
%    N.B. The diffraction-limited model is a special case of the
%    shift-invariant model with the PSF constrained to be the ideal blur
%    function of a circular aperture.
%
%    About the OTF and PSF
%       We store the OTF data with DC at (1, 1). This is true throughout
%       isetbio. To understand the implications for certain calculations
%       see t_codeFFTinMatlab.
%
%    Although we use the Matlab-style representation (DC at (1, 1)) for the
%    OTF, when we make graphs and images we put the center of the image at
%    the center -- of course -- and we also put the DC value of the OTF in
%    the middle of the image. Hence, when we return the frequency support
%    or the spatial support we create values for frequencies that run from
%    negative to positive with 0 sf in the middle. Similarly, when we
%    compute the spatial support we create spatial samples that run below
%    and above zero. If we've done things correctly, 0 sf should at
%    location floor(N / 2) + 1, where N is the number of frequency samples.
%    To convert the OTF so that it matches up with this representation,
%    apply fftshift to it.
%
%    To get the PSF from the OTF, we use the PTB routine OtfToPsf. This
%    expects the DC term in the middle, so we call
%
%       [~, ~, PSF] = OtfToPsf([], [], fftshift(OTF))
%
%    to do the converstion. The null args to OtfToPsf involve the support
%    and are not needed here.
%
% Inputs:
%    optics   - Struct. An optics structure from which to retrieve the
%               parameter value.
%    parm     - String. The parameter you wish to retrive (or calculate)
%               the value of. There are a large number of options, so they
%               are broken into categories below. Note, an asterisk (*)
%               next to the parameter name indicates that you can specify
%               the parameter units using varargin. When available, the
%               type of the value is specified.
%        Bookkeeping
%           {'name'}                - String. Name for these optics
%           {'type'}                - String. Always 'optics'
%           {'model'}               - String. Type of optics computation,
%                                     diffractionLimited, rayTrace, or
%                                     shiftInvariant.
%           {'fnumber'}             - Numeric. f#, ratio of focal length to
%                                     aperture, a dimensionless quantity.
%           {'effective fnumber'}   - Numeric. Effective f-number
%           {'focal length'}        - Numeric. Focal length (M)
%           {'power'}               - Numeric. Optical power in diopters
%                                     (1/f), units 1/M
%           {'image distance'}      - Numeric. Image distance from
%                                     lensmaker's equation.
%           {'image height'}*       - Numeric. Image height
%           {'image width'}*        - Numeric. Image width
%           {'image diagonal'}*     - Numeric. Image diagonal size
%           {'numerical aperture'}  - Numeric. Numerical aperture
%           {'aperture diameter'}*  - Numeric. Aperture diameter
%           {'aperture radius'}*    - Numeric. Aperture radius
%           {'aperture area'}*      - Numeric. Aperture area
%           {'magnification'}       - Numeric. Optical image magnification
%                                     (<0 inverted)
%           {'pupil magnification'} - Numeric. Magnification for pupil.
%        Off-axis methods and data
%           {'off axis method'} - String. The custom relative illumination
%                                 method to use.
%           {'cos4th method'}   - default cos4th method
%           {'cos4th data'}     - Matrix. place to store cos4th data
%           {'vignetting'}      - Deprecated.  Left from old sceneEye code
%
%        OTF information - Used for shift-invariant calculations
%           {'otf data'}     - Matrix. the optical transfer function data
%           {'otf size'}     - Vector. 1x2 of x and y sizes.
%           {'otf fx'}       - Vector. column (fx) samples of OTF data
%           {'otf fy'}       - Vector. row (fy) samples of OTF data
%           {'otf support'}  - Cell array. val{1:2}, of fy, fx samples
%           {'otf wave'}     - Vector. wavelength samples of the otf data
%           {'otf binwidth'} - Numeric. difference in wavelength samples
%           {'psf data'}     - Matrix. psf data, which is calculated from
%                              the stored otfdata
%           {'diffraction limited psf data'}*
%                            - Matrix. diffraction limited psf data,
%                              requires an additional argument of
%                              wavelength, with optional arguments for
%                              units, the number of samples, and the
%                              frequency over the samples.
%           {'psf spacing'}  - Numeric. PSF Spacing.
%           {'psf support'}  - Array. 1x2 Cell array.
%           {'incoherentcutoffspatialfrequency'} *
%                            - Vector. A Vector of incoherent cutoff freq
%                              for all of the wavelengths
%           {'maxincoherentcutoffspatialfrequency'} *
%                            - Numeric. Largest incoherent cutoff
%        Wavlength information
%           {'spectrum'}      - wavelength information
%               {'wavelength'} - Vector. wavelength samples
%               {'nwave'}      - Numeric. number of wavelength samples
%               {'binwidth'}   - Numeric. spacing between the samples
%
%           For the diffraction limited case, we store transmittance
%           in these slots.
%           {'transmittance'} - Vector. Transmittance function of the lens
%               {'wave'}       - Wavelength samples *lenstransmittancewave*
%               {'scale'}      - Vector. Spectral radiance scale factor
%                                *lensTransmittanceScale*
%           For human optics, we store a human lens object and derive
%           the lens transmittance and other properties on the fly
%           from human data. Only the lens pigment density can be
%           changed.
%           {'lens'}          - Struct. The lens object
%
%    varargin - (Optional). Additional arguments for retrieving or
%               calculating the desired parameter. Usually the units for
%               the calculations.
%
% Outputs:
%    val      - VARIES. The value of the requested parameter.
%
% Optional key/value pairs:
%     None.
%
%
% Notes:
%    * [Note - DHB:  See notes in code below about usage for
%      'diffractionlimitedpsfdata' and corresponding frequency support
%      calls. These are only used in oiPlot to plot a diffraction limited
%      function, and should be used with caution. We may want to re-write
%      to be more consistent in our usage across different ways of getting
%      the dl information, but this would be fairly involved and may not be
%      worth it.]
%    * [Note - DHB: Not all options appear to be documented in the above
%      header comments.]
%    * TODO: When 2015B phased out, replace strfind() with contains()
%    * Data representation: the OTF and PSF
%         We store the OTF data with DC at (1, 1). This is true throughout
%         isetbio. To understand the implications for certain calculations
%         see t_codeFFTinMatlab.
%
%         Although we use the Matlab-style representation of (DC at (1, 1))
%         for the OTF, when we make graphs and images we put the center of
%         the image at the center (of course) and we also put the DC value
%         of the OTF in the middle of the image. Hence, when we return the
%         frequency support or the spatial support we create values for
%         frequencies that run from negative to positive with 0 sf in the
%         middle. Similarly, when we compute the spatial support we create
%         spatial samples that run below and above zero. If we've done
%         things correctly, 0 sf should be at location floor(N/2)+1, where
%         N is the number of frequency samples. To convert the OTF so that
%         it matches up with this representation, apply an fftshift to it.
%
%         To get the PSF from the OTF, we use the PTB routine OtfToPsf.
%         This expects the DC term in the middle, so we call
%
%               [~, ~, PSF] = OtfToPsf([], [], fftshift(OTF))
%
%         to do the conversion. The null args to OtfToPsf involve the
%         support and are not needed here.
%
% Copyright Imageval, LLC, 2005
%
% See Also:
%   oiCreate, oiGet, opticsSet
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    12/21/17  dhb  Use OtfToPsf to do the conversion, but with backwards
%                   compatible control.
%    03/14/18  jnm  Formatting
%    08/30/18  dhb  Replace abs(fft2(OTF)) with OtfToPsf for M-W optics,
%                   which was special cased for backward compatibility.
%                   This is based on the view that what was being done was
%                   not correct, and that we should bite the bullet and
%                   change. See issue #373 on github.
%    10/18/18  dhb  Special case OtfToPsf negative tolerance for M-W optics
%                   versus everything else. Be explicit about tolerance in
%                   all the calls, to remove sensitivity to OtfToPsf
%                   default. This fixes a broken example.
%              dhb  Deleted a pile of commented out code.
%    06/26/19  JNM  Documentation update & some formatting fixes.

% Examples:
%{
    oi = oiCreate;
    oi = oiSet(oi, 'wave', 400:10:700);
    NA = oiGet(oi, 'optics na');              % Numerical aperture
    psf = oiGet(oi, 'optics psf data', 600);  % Shift invariant data
    sSupport = oiGet(oi, 'optics psf support');
    vcNewGraphWin;
    mesh(sSupport{1}, sSupport{2}, psf);
    otf = oiGet(oi, 'optics otf data', 450);
    vcNewGraphWin;
    mesh(fftshift(abs(otf)));
%}
%{
    oi = oiCreate;
    % The direct calls using opticsGet are:
    optics = oiGet(oi, 'optics');
    otf450 = opticsGet(optics, 'otf data', 450);
    otfSupport = opticsGet(optics, 'otf support');  % Cycles/mm
    vcNewGraphWin;
    mesh(otfSupport{1}, otfSupport{2}, fftshift(abs(otf450)))
%}
%{
    % Check that dist per deg computation is correct.
    oi = oiCreate('human');
    distPerDeg = oiGet(oi,'optics dist per deg')*10^3; % mm/deg
    degPerDist = oiGet(oi,'optics deg per dist')*10^-3; % deg/mm
    % Check these values
    assert(abs(distPerDeg-0.3)<0.1); % close to 0.3 mm/deg
    assert(abs(degPerDist-3.3)<0.1); % close to 3.3 deg/mm
%}

%% Control some printout
RESPECT_THE_COMMAND_WINDOW = true;

val = [];
if ~exist('optics', 'var') || isempty(optics)
    error('No optics specified.');
end
if ~exist('parm', 'var') || isempty(parm)
    error('No parameter specified.');
end

parm = ieParamFormat(parm);
switch parm
    case 'name'
        val = optics.name;

    case 'type'
        val = optics.type;  % Should always be 'optics'

    case {'fnumber', 'f#'}
        % This is the f# assuming an object is infinitely far away.
        val = optics.fNumber;

    case {'model', 'opticsmodel'}
        if checkfields(optics, 'model'), val = optics.model;
        else, val = 'diffraction limited';
        end

        % The user can set 'Diffraction limited' but have
        % 'diffractionlimited' returned.
        val = ieParamFormat(val);

    case {'effectivefnumber', 'efffnumber', 'efff#'}
        % The f# if the object is not at infinity.
        val = opticsGet(optics, 'fNumber') * ...
            (1 - opticsGet(optics, 'mag'));

    case {'focallength', 'flength'}
        % opticsGet(optics, 'flength', units);
        val = optics.focalLength;
        if ~isempty(varargin)
            val = val * ieUnitScaleFactor(varargin{1});
        end

    case {'power', 'diopters'}
        % opticsGet(optics, 'power', 'mm')
        % Diopters (1/m).
        % Sometimes we ask for a different unit when we don't want
        % diopters, but we want dist per deg in another unit.
        units = 'm';
        if ~isempty(varargin), units = varargin{1}; end
        val = 1 / opticsGet(optics, 'focallength', units);

    case {'imagedistance', 'focalplane', 'focalplanedistance'}
        % Lensmaker's equation calculation of focal plane distance from
        % the center of the lens. If no source distance is provided, we
        % assume infinite distance
        %
        % opticsGet(optics, 'focalplane', sDist);  % sDist = sourceDistance
        % opticsGet(optics, 'focalplanedistance'); % Infinite source dist

        fL = opticsGet(optics, 'focal Length');
        if isempty(varargin), sDist = Inf;
        else, sDist = varargin{1};
        end
        val = 1 / (1/fL - 1/sDist);   % Lens maker equation

    % These all need the scene field of view. They can be computed from the
    % geometry of the image distance and FOV.
    case {'imageheight'}
        % opticsGet(optics, 'imageheight', fov) -- field of view in degrees
        if isempty(varargin)
            disp('fov required.');
            return;
        else
            fov = varargin{1};
            imageDistance = opticsGet(optics, 'focal plane distance');
            val = 2 * imageDistance * tand(fov / 2);
            if length(varargin) < 2, return;
            else, val = ieUnitScaleFactor(varargin{2}) * val;
            end
        end

    case {'imagewidth'}
        % opticsGet(optics, 'imagewidth', fov) -- fov in degrees
        % opticsGet(optics, 'imagewidth', fov, 'mm') -- fov in degrees,
        % output units in mm
        if isempty(varargin), return;
        else
            fov = varargin{1};
            imageDistance = opticsGet(optics, 'focalplanedistance');
            val = 2 * imageDistance * tand(fov / 2);
            if length(varargin) < 2, return;
            else, val = ieUnitScaleFactor(varargin{2}) * val;
            end
        end

    case {'imagediagonal', 'diagonal'}
        % opticsGet(optics, 'imagediagonal', fov)  -- fov in degrees
        if isempty(varargin), return;
        else
            fov = varargin{1};
            h = opticsGet(optics, 'imageheight', fov, 'm');
            w = opticsGet(optics, 'imagewidth', fov, 'm');
            val = sqrt(h ^ 2 + w ^ 2);
        end
        if length(varargin) < 2, return;
        else, val = ieUnitScaleFactor(varargin{2}) * val;
        end

    case {'na', 'numericalaperture'}
        % Should this be a call to effective f#?
        val = 1 / (2 * opticsGet(optics, 'fnumber'));

    case {'aperturediameter', 'diameter', 'pupildiameter'}
        val = opticsGet(optics, 'focalLength') / ...
            opticsGet(optics, 'fnumber');
        if ~isempty(varargin)
            val = ieUnitScaleFactor(varargin{1}) * val;
        end

    case {'apertureradius', 'radius', 'pupilradius'}
        val = opticsGet(optics, 'diameter')/2;
        if ~isempty(varargin)
            val = ieUnitScaleFactor(varargin{1}) * val;
        end

    case {'aperturearea', 'pupilarea'}
        val = pi * ((opticsGet(optics, 'focalLength') / ...
            opticsGet(optics, 'fnumber')) * (1 / 2)) ^ 2;
        if ~isempty(varargin)
            val = ieUnitScaleFactor(varargin{1}) ^ 2 * val;
        end

    case {'magnification', 'mag'}
        % If the ray trace magnification is present, use that. Otherwise
        % make a magnification estimate using the source distance and focal
        % length (via lensmaker equation).
        % opticsGet(optics, 'mag', sDist) -- specify source distance
        % opticsGet(optics, 'mag')        -- current source distance

        % Skip model has unit magnification, but still negative
        if strcmpi(opticsGet(optics, 'model'), 'skip')
            val = -1;
            return;
        elseif length(varargin) < 1
            scene = vcGetObject('scene');
            if ~isempty(scene)
                sDist = sceneGet(scene, 'distance');
            else
                % What should we do here? No scene and no distance
                % specified. Do we make the magnification 0? Seems OK.
                val = 0;
                return;
            end
        else
            sDist = varargin{1};
        end
        % Distance to the image divided by distance to the object (source)
        % http://en.wikipedia.org/wiki/Magnification
        % Also f / (f - distObj)
        % If you want the mag to be -2, then
        % distObj -2 * f + 2 * distObj
        %   = f, 2 * distObj
        %   = 3 * f, distObj = 3f/2
        % and in general distObj = (-M+1)f/(-M) * f
        val = -(opticsGet(optics, 'focalPlaneDistance', sDist))/sDist;

    case {'pupilmagnification', 'pupmag'}
        % Pupil magnification is the ratio of the exit pupil to the input
        % pupil (diameters)
        val = 1;

    case {'wave', 'wavelength', 'wavelengthvalues'}
        % Wavelength related
        % nanometers is default
        % opticsGet(optics, 'wavelength', unit)
        %
        if checkfields(optics, 'spectrum', 'wave')
            val = optics.spectrum.wave;
        end
        if isempty(val)
            scene = vcGetObject('scene'); val = sceneGet(scene, 'wave');
            if ~isempty(val)
                if (~RESPECT_THE_COMMAND_WINDOW)
                    disp([' * * Optics wave sampling set by selected ' ...
                        'scene.']);
                end
            end
        end
        if isempty(val)
            val = 400:10:700;
            val = val(:);
            if (~RESPECT_THE_COMMAND_WINDOW)
                disp(' * * Optics wave sampling set by selected scene.');
            end
        end

        if ~isempty(varargin)
            s = ieUnitScaleFactor(varargin{1}) / ieUnitScaleFactor('nm');
            val = val * s;
        end

    case {'nwave', 'numberofwavelengthsamples'}
        if checkfields(optics, 'spectrum', 'wave')
            val = length(optics.spectrum.wave);
        end

    case {'binwidth', 'wavelengthbinwidth'}
        % Nanometers
        wave = opticsGet(optics, 'wave');
        if length(wave) > 1, val = wave(2) - wave(1);
        else, val = 1;
        end

    case {'transmittance', 'lenstransmittance', 'lenstransmittancescale'}
        % Lens spectral transmittance representation
        % opticsGet(optics, 'transmittance', wave)
        %
        % When the optics are human, there is a lens object and we read the
        % transmittance from it. When the optics are diffraction, there is
        % a transmittance.wave and transmittance.scale field, and we
        % interpolate from that.

        if checkfields(optics, 'lens')  % Human optics
            % The lens gets the proper wavelength sample because the object
            % has a 'listener'
            val = optics.lens.get('transmittance', varargin{:});

        elseif checkfields(optics, 'transmittance') % Diffraction optics
            val = optics.transmittance.scale;

            % No listener for the old-fashioned struct.
            if ~isempty(varargin)
                newWave = varargin{1};
                wave = optics.transmittance.wave;
                scale = optics.transmittance.scale;

                if min(newWave(:))< min(wave(:)) || ...
                        max(newWave(:)) > max(wave(:))
                    % Extrapolation required.
                    disp('Extrapolating transmittance with 1''s')
                    val = interp1(wave, scale, newWave, 'linear', 1)';
                else
                    val = interp1(wave, scale, newWave, 'linear')';
                end
            end
        end

    case {'transmittancewave', 'lenstransmittancewave'}
        % The lens transmittance wavelength samples
        % Different case for lens and for
        if checkfields(optics, 'lens')  % Human optics
            % The lens gets the proper wavelength sample because the object
            % has a 'listener'
            val = optics.lens.get('wave');
        elseif checkfields(optics, 'transmittance')
            val = optics.transmittance.wave;
        else
            error('No lens or transmittance');
        end

    case {'lens'}
        % New lens object.
        val = optics.lens;

        % ----- Diffraction limited parameters
    case {'dlfsupport', 'dlfsupportmatrix'}
        % Two different return formats. Either
        %  val{1} and val{2} as vectors, or
        %  val = fSupport(:, :, :);
        % opticsGet(optics, 'dl fsupport', wave, unit, nSamp)
        % opticsGet(optics, 'dl fsupport matrix', wave, unit, nSamp)
        %
        % Diffraction limited frequency support at a wavelength (i.e.,
        % support out to the incoherent cutoff frequency). This can be
        % used for plotting, for example
        %
        % * [Note - DHB: The returned array has dimension = 2 * nSamp. This
        %   is confusing, to me at least.]
        if length(varargin) < 1, error('Must specify wavelength');
        else, thisWave = varargin{1}; end
        if length(varargin) < 2, units = 'mm';
        else, units = varargin{2}; end
        if length(varargin) < 3, nSamp = 30; else, nSamp = varargin{3}; end

        % Sometimes the optics wavelength hasn't been defined because, say,
        % we haven't run through a scene. So we trap that case here.
        waveList = opticsGet(optics, 'wavelength');
        idx = ieFindWaveIndex(waveList, thisWave);
        inCutoff = opticsGet(optics, 'inCutoff', units);
        inCutoff = inCutoff(idx);

        fSamp = (-nSamp:(nSamp - 1)) / nSamp;
        val{1} = fSamp * inCutoff;
        val{2} = fSamp * inCutoff;

        % Alternative return format
        if strfind(parm, 'matrix')
            [valMatrix(:, :, 1), valMatrix(:, :, 2)] = ...
                meshgrid(val{1}, val{2});
            val = valMatrix;
        end

    case {'incoherentcutoffspatialfrequency', 'incutfreq', 'incutoff'}
        % cycles/distance
        % Cutoff spatial frequency for a diffraction limited lens. See
        % formulae in dlCore.m
        apertureDiameter = opticsGet(optics, 'aperture diameter');
        imageDistance = opticsGet(optics, 'focal plane distance');
        wavelength = opticsGet(optics, 'wavelength', 'meters');

        % Sometimes the optics wavelength have not been assigned because
        % there is no scene and no oiCompute has been run. So, we can just
        % choose a sample set.
        if isempty(wavelength), wavelength = (400:10:700) * 10 ^ -9; end

        % See dlCore.m for a description of the formula. We divide by the
        % scale factor, instead of multiplying, because these are
        % frequencies (1/m), not distances.
        val = (apertureDiameter / imageDistance) ./ wavelength;
        if ~isempty(varargin)
            val = val / ieUnitScaleFactor(varargin{1});
        end

    case {'maxincoherentcutoffspatialfrequency', ...
            'maxincutfreq', 'maxincutoff'}
        % opticsGet(optics, 'maxincutoff', 'm')
        % opticsGet(optics, 'maxincutoff')
        if isempty(varargin), val = max(opticsGet(optics, 'incutoff'));
        else, val = max(opticsGet(optics, 'incutoff', varargin{1}));
        end

        % -------   OTF information and specifications. This case is used
        % for shift-invariant calculations. The ray trace structures below
        % are used for non shift-invariant cases derived from Zemax or Code
        % V data sets.

    case {'otf', 'otfdata', 'opticaltransferfunction'}
        % You can ask for a particular wavelength with the syntax
        %
        %    opticsGet(optics, 'otfData', oi, spatialUnits, wave)
        %
        % OTF values can be complex. They are related to the PSF data by
        % the fft. We use PstToOtf to do this, so we use consistent
        % conventions throughout isetbio. But, we still have to be careful
        % about whether the OTF is zero centered or centered with DC term
        % at the upper left.
        opticsModel = opticsGet(optics, 'model');
        thisWave = [];
        switch lower(opticsModel)
            case 'diffractionlimited'
                % For diffraction limited case, the call must be
                % otf = opticsGet(...
                %    optics, 'otf data', oi, fSupport, [wave]);
                units = 'mm';
                if isempty(varargin)
                    error('format is ... oi, spatialUnits, wave');
                else, oi = varargin{1};
                end

                if length(varargin) > 1, units = varargin{2}; end
                if length(varargin) > 2, thisWave = varargin{3}; end

                fSupport = oiGet(oi, 'fSupport', units);  % 'cycles/mm'
                wavelength = oiGet(oi, 'wave');
                OTF = dlMTF(oi, fSupport, wavelength, 'millimeters');

            case 'shiftinvariant'
                % opticsGet(optics, 'otf data', [wave]);
                if checkfields(optics, 'OTF', 'OTF')
                    OTF = optics.OTF.OTF;
                    if ~isempty(varargin), thisWave = varargin{1}; end
                else, OTF = [];
                end

            case 'raytrace'
                error(['opticsGet(optics, ''OTF'') not supported for ' ...
                    'ray trace']);

            otherwise
                error('OTFData not implemented for %s model', opticsModel);
        end

        % Wavelength is asked for
        if ~isempty(thisWave)
            [idx1, idx2] = ...
                ieWave2Index(opticsGet(optics, 'otfWave'), thisWave);
            if idx1 == idx2
                val = OTF(:, :, idx1);
            else
                wave = opticsGet(optics, 'otfwave');
                w = 1 - ((varargin{1} - wave(idx1)) / ...
                    (wave(idx2) - wave(idx1)));
                val = (w * OTF(:, :, idx1) + (1 - w) * OTF(:, :, idx2));
                % wave(idx1), varargin{1}, wave(idx2), w)
            end
        else
            % Returns the entire data set
            val = OTF;
        end

    case {'degreesperdistance', 'degperdist'}
        % opticsGet(optics, 'deg per dist', 'mm')
        % We use this constant to convert from the input spatial frequency
        % units (cycles/deg) to cycles/meter needed for the Hopkins eye. We
        % need to calculate this value from the optics data passed in.
        %
        % Given D0, the focal plane is 1/D0 meters from the lens. Call this
        % the adjacent edge of the right triangle from the image plane to
        % the lens.
        %
        % 1 deg of visual angle is
        %   atan(opp / (1 / D0)) = 1             (deg)
        %   opp / (1 / D0) = tand(1)           (1 / rad)
        %   opp = tand(1) * (1 / D0)           (1 / rad * meter)
        %   1 / opp = 1 / (tand(1) * (1 / D0)) (rad / meter)
        %
        % The conversion: (cycles / rad) * (rad / meter) = cycles / meter
        units = 'm';
        if ~isempty(varargin), units = varargin{1}; end
        D0 = opticsGet(optics, 'power', units);
        val = 1 /(1 / D0 * tand(1));

    case {'distperdeg', 'distanceperdegree'}
        units = 'm';
        if ~isempty(varargin), units = varargin{1}; end
        val = 1 / opticsGet(optics, 'deg per dist', units);

    %------------------------
    % [Note: DHB - I think we should delete the stuff below and the OTF.fx
    % entry, and only use oiFrequencySupport. This will require debugging,
    % particularly in the script s_TestSI and perhaps other places.]
    case {'otfsupport'}
        % val = opticsGet(optics, 'otf support', 'mm');
        %
        % Row & col (Y, X) spatial frequency range for the OTF data [Y, X]
        % meshgrid(val{1}, val{2}) produces the matrices for surface
        % plotting. Frequency is stored in non-standard units of cycles/mm.
        % This will be annoying to fix some day, sigh.
        units = 'mm';
        if ~isempty(varargin), units = varargin{1}; end
        val{1} = opticsGet(optics, 'otf fy', units);
        val{2} = opticsGet(optics, 'otf fx', units);

    case {'otffx'}
        % cycles/mm!!! Non-standard unit. Must fix up some day.
        if checkfields(optics, 'OTF', 'fx'), val = optics.OTF.fx(:)'; end
        % Transform into other units if required
        if ~isempty(varargin)
            unit = ieParamFormat(varargin{1});
            if strcmp(unit, 'cycles/deg') || strcmp(unit, 'cyclesperdeg')
                val = val * tand(1) * ...
                    opticsGet(optics, 'focal length', 'mm');
            else
                val = (val * 10 ^ 3) / ieUnitScaleFactor(unit);
            end
        end

    case {'otffy'}
        % cycles/mm!!! Non-standard unit. TODO: Must fix up some day.
        if checkfields(optics, 'OTF', 'fy'), val= optics.OTF.fy(:)'; end
        % Put into meters and then apply scale factor
        if ~isempty(varargin)
            unit = ieParamFormat(varargin{1});
            if strcmp(unit, 'cycles/deg') || strcmp(unit, 'cyclesperdeg')
                val = val * tand(1) * ...
                    opticsGet(optics, 'focal length', 'mm');
            else
                val = (val * 10 ^ 3) / ieUnitScaleFactor(unit);
            end
        end

    case {'otfsize'}
        % Row and col samples
        if checkfields(optics, 'OTF', 'OTF')
            tmp = size(optics.OTF.OTF); val = tmp(1:2);
        end

    case {'otfwave'}
        % opticsGet(optics, 'otf wave', 'nm');
        % nm is the default.
        % This should probably go away and we should only wave 'wave'.
        if checkfields(optics, 'OTF', 'wave'), val = optics.OTF.wave;
        else, val = opticsGet(optics, 'wave');

        end
        if ~isempty(varargin)
            units = varargin{1};
            val = val * 1e-9 * ieUnitScaleFactor(units);
        end

    case {'otfbinwidth'}
        otfWave = opticsGet(optics, 'otfWave');
        if length(otfWave)>1, val = otfWave(2) - otfWave(1);
        else, val = 1;
        end

    case {'psfdata'}
        % To return the psf at 500 nm use
        %
        %    psf = opticsGet(optics, 'psfData', 500);
        %    mesh(psf);
        %
        % The diffraction limited and shift invariant models are
        % handled slightly differently. The dlMTF is needed for the
        % diffraction case, while the OTF holds the data in the
        % shift-invariant case.
        if length(varargin) < 1, thisWave = opticsGet(optics, 'wave');
        else, thisWave = varargin{1};
        end
        nWave = length(thisWave);

        if length(varargin) < 2, units = 'um';
        else, units = varargin{2};
        end

        oModel = opticsGet(optics, 'model');
        switch lower(oModel)
            case 'diffractionlimited'
                % opticsGet(optics, 'psf Data', 500, 'um');
                % nSamp - # of frequency steps from 0 to incoherent cutoff
                nSamp = 100;
                dlF = opticsGet(optics, 'dlFSupport', thisWave(1), ...
                    units, nSamp);
                [fSupport(:, :, 1), fSupport(:, :, 2)] = ...
                    meshgrid(dlF{1}, dlF{2});

                % Get diffraction limited OTF and derive PSF from it.
                otf = dlMTF(optics, fSupport, thisWave, units);
                if nWave == 1
                    [~, ~, val] = OtfToPsf([], [], fftshift(otf), ...
                        'negativeFractionalTolerance', 1e-3);
                else
                    val = zeros(size(otf, 1), size(otf, 2), ...
                        length(thisWave));
                    for ii = 1:length(thisWave)
                        [~, ~, val(:, :, ii)] = OtfToPsf([], [], ...
                            fftshift(otf(:, :, ii)), ...
                            'negativeFractionalTolerance', 1e-3);
                    end
                end

            case 'shiftinvariant'
                if checkfields(optics, 'OTF', 'OTF')
                    if nWave == 1
                        % Just do one specified wavelength
                        otf = opticsGet(optics, 'otf data', thisWave);
                        if (strcmp(opticsGet(optics, 'name'), 'human-MW'))
                            [~, ~, val] = OtfToPsf([], [], ...
                                fftshift(otf), ...
                                'negativeFractionalTolerance', 5e-2);
                        else
                            [~, ~, val] = OtfToPsf([], [], ...
                                fftshift(otf), ...
                                'negativeFractionalTolerance', 1e-3);
                        end
                    else
                        % Do all the wavelenghts
                        val = zeros(size(optics.OTF.OTF));
                        for ii = 1:length(thisWave)
                            if strcmp(opticsGet(optics, 'name'), ...
                                    'human-MW')
                                [~, ~, val(:,:,ii)] = OtfToPsf([], [], ...
                                    fftshift(optics.OTF.OTF(:, :, ii)), ...
                                    'negativeFractionalTolerance', 5e-2);
                            else
                                [~, ~, val(:,:,ii)] = OtfToPsf([], [], ...
                                    fftshift(optics.OTF.OTF(:, :, ii)), ...
                                    'negativeFractionalTolerance', 5e-2);
                            end
                        end
                    end

                    % PSF should be real, make it so but warn if it is not.
                    if ~isreal(val)
                        warning('ISET:complexpsf', 'complex psf');
                        val = abs(val);
                    end
                else
                    warning(['ISET:otfdata', 'No OTF data stored in ' ...
                        'optics.'])
                end
        end

    case {'psfspacing'}
        % opticsGet(optics, 'psf spacing', unit)
        % Sample spacing of the psf points
        %
        % Warning:  We assume that fx and fy have the same peak spatial
        % frequency and spatial sampling.

        if ~isempty(varargin), units = varargin{1}; else, units = 'mm'; end

        % We use these parameters for the diffraction limited case.
        if length(varargin) > 1, thisWave = varargin{2}; end
        if length(varargin) > 2, nSamp = varargin{3}; end

        oModel = opticsGet(optics, 'model');
        switch oModel
            case 'shiftinvariant'
                fx = opticsGet(optics, 'otf fx', units);
                if isempty(fx)
                    error('No otffx calculated yet. Fix me.');
                end
            case 'diffractionlimited'
                fx = opticsGet(optics, 'dl fSupport', thisWave(1), ...
                    units, nSamp);
                fx = fx{1};
            otherwise
                error('Unknown optics model %s\n', oModel);
        end

        % Peak frequency in cycles/distance (units). We have two samples in
        % that distance, so the sample spacing is half that distance.
        peakF = max(fx(:));
        val = 1 / (2 * peakF);

    case {'psfsupport'}
        % opticsGet(optics, 'psf support', unit)
        % Calculate the spatial support for the psf.
        %
        % Returns mesh grid of X & Y values in a 2-dimensional cell array.
        % Used for mesh plotting often. X/Y could be reversed in 1 and 2.

        % These should be integrated and coordinated with fSupport. If we
        % would like the spatial support to be smaller and finer, we should
        % scale fSupport * 4.
        %{
           samp = (-nSamp:(nSamp - 1));
           [X, Y] = meshgrid(samp, samp);
           deltaSpace = 1 / (2 * max(fSupport(:)));
           sSupport(:, :, 1) = X * deltaSpace;
           sSupport(:, :, 2) = Y * deltaSpace;
        %}

        if length(varargin) >= 1, units = varargin{1};
        else, units = 'mm'; end
        oModel = opticsGet(optics, 'model');

        switch oModel
            case 'shiftinvariant'
                sz = opticsGet(optics, 'otf size');
                if isempty(sz), error('No optical image data'); end
                if (sz(1) ~= sz(2)), error('OTF support not square'); end

                % Create one-dimensional integer samples.
                % Handle case of sz even versus sz odd. I am more confident
                % of the case where n is even.
                if rem(sz(1), 2) == 0
                    n = sz(1) / 2;
                    x = -n:(n - 1);
                else
                    n = floor(sz(1) / 2);
                    x = -n:n;
                end
                x = x * opticsGet(optics, 'psf spacing', units);

            case 'diffractionlimited'
                wave = opticsGet(optics, 'wave');
                thisWave = wave(1);
                psf  = opticsGet(optics, 'psf data', thisWave);
                if isempty(psf), error('No psf data'); end
                sz = size(psf);
                if (sz(1) ~= sz(2)), error('psf support not square'); end

                % Create one-dimensional integer samples. Handle case of sz
                % even versus sz odd. I am more confident of the case
                % where n is even.
                if rem(sz(1), 2) == 0
                    n = sz(1) / 2;
                    x = -n:(n - 1);
                else
                    n = floor(sz(1) / 2);
                    x = -n:n;
                end
                nSamp = length(x);
                x = x * opticsGet(optics, 'psf spacing', units, ...
                    thisWave, nSamp);

            otherwise
                error('unknown model %s\n', oModel)
        end

        % Make grids
        [X, Y] = meshgrid(x, x);
        val{1} = X;
        val{2} = Y;

    %----------- Relative illumination (off-axis) specifications
    case {'offaxis', 'offaxismethod', 'relativeilluminationtype'}
        % This is the method used to compute relative illumination. It can
        % be 'Skip', 'cos4th', 'codeV', or 'Zemax'.
        val = optics.offaxis;

    case {'cos4thmethod', 'cos4thfunction'}
        % Most people use cos4th as an offaxis method. In that case, the
        % function that implements cos4th can be stored here. I suspect
        % this extra step is not needed. We have a cos4th function and we
        % use that without allowing some other implementation. It is here
        % only for some old backwards compatibility.
        %
        % Do not run cos4th when you are using the Code V (or probably
        % Zemax methods. These calculations include the cos4th mechanisms
        % in the lens calculations.

        if checkfields(optics, 'cos4th', 'function')
            val = optics.cos4th.function;
        end
    case {'cos4th', 'cos4thdata', 'cos4thvalue'}
        % Numerical values. Should change field to data from value. I don't
        % think this is ever used, is it?
        if checkfields(optics, 'cos4th', 'value')
            val = optics.cos4th.value;
        end

    otherwise
        error('Unknown optics parameter %s\n', parm);
end

end