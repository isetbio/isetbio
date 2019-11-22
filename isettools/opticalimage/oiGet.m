function val = oiGet(oi, parm, varargin)
% Get properties and derived quantities from an optical image structure
%
% Syntax:
%   val = oiGet(oi, parm, [varargin])
%
% Description:
%    Optical image parameters are stored in a structure. Some of the
%    parameters are stored directly, and others are calculated from the
%    data structure. We store the unique values and calculate many derived
%    values in this routine.
%
%    The optical image structures are stored in a cell array,
%    vcSESSION.OPTICALIMAGE{}
%
%    To retrieve the currently selected optical image, use either
%        oi = vcGetObject('OI');
%        [val, oi] = vcGetSelectedObject('OI');
%
%    A '*' indicates that the syntax sceneGet(scene, param, unit) can be
%    used, where unit specifies the spatial scale of the returned value:
%    'm', 'cm', 'mm', 'um', 'nm'. Default is meters ('m').
%
%    The optics data structure is stored in the oi as oi.optics. It is
%    possible to retrieve the optics parameters as
%
%        v = oiGet(oi, 'optics fnumber');
%
%    The lens data structure is stored as oi.optics.lens. You can query the
%    lens properties using
%
%        v = oiGet(oi, 'lens property');
%
%     To see the full range of optics parameters, use doc opticsGet.
%
%    There are examples contained in the code. To access the examples, type
%    'edit oiGet.m' into the Command Window.
%
% Inputs:
%    oi       - Struct. An optical image structure.
%    parm     - String. The parameter you wish to retrieve the info from.
%               Some of the options (in categories) include the following:
%          General properties
%               {'name'}            - String. The optical image name.
%               {'type'}            - String. Always 'opticalimage'.
%               {'filename'}        - String. If read from file, put here.
%               {'consistency'}     - Boolean. Whether the oiWindow display
%                                     reflects the current state 1 (yes) or
%                                     0 (not).
%               {'rows'}            - Numeric. The number of row samples.
%               {'cols'}            - Numeric. The number of col samples.
%               {'size'}            - Matrix. [rows, cols].
%               {'image distance'}  - Numeric. The distance from lens to
%                                     image, negative.
%               {'hfov'}            - Numeric. The horizontal field of view
%                                     in degrees.
%               {'vfov'}            - Numeric. The vertical field of view
%                                     in degrees.
%               {'aspectratio'}     - Numeric. The aspect ratio of image.
%               {'height'}*         - Numeric. The image height.
%               {'width'}*          - Numeric. The image width.
%               {'diagonal'}*       - Numeric. The image diagonal length.
%               {'heightandwidth'}* - Matrix. [height, width].
%               {'area'}*           - Numeric. The optical image area.
%               {'centerpixel'}     - Matrix. The [row, col] point at the
%                                     image center.
%          Irradiance
%               {'data'}                - Struct. The data structure.
%                  {'photons'}          - Matrix. Irradiance data
%                  {'photons noise'}    - Matrix. Irradiance data w/ photon
%                                         noise for 50 msec integration
%                                         time & 1 pixel collection area.
%                  {'data max'}         - Numeric. Used for compression,
%                                         not for general users.
%                  {'data min''}        - Numeric. Used for compression,
%                                         not for general users.
%                  {'bit depth''}       - Numeric. Used for compression,
%                                         not for general users.
%                  {'energy'}           - Matrix. Energy rather than
%                                         photon representation.
%                  {'roi energy'}       - Numeric. Energy of the points in
%                                         a region of interest
%                  {'roi mean energy'}  - Numeric. Average energy of the
%                                         points in a region of interest.
%                  {'energy noise'}     - Matrix. Energy with photon noise
%                                         for 50 msec integration time and
%                                         1 pixel collection area.
%                  {'mean illuminance'} - Numeric. The mean illuminance.
%                  {'illuminance'}      - Matrix. A matrix representation
%                                         of the optical image illuminance.
%                  {'xyz'}              - Matrix. [row, col, 3] an image of
%                                         the irradiance XYZ values.
%                  {'lms'}              - Matrix. A 3D array of cone
%                                         values, from Stockman.
%          Wavelength information
%               {'spectrum'}    - Struct. Wavelength information structure.
%                  {'binwidth'} - Numeric. The spacing between samples.
%                  {'wave'}     - Array. Wavelength samples (nm).
%                  {'nwave'}    - Numeric. Number of wavelength samples.
%          Resolution
%               {'hspatial resolution'}*  - Numeric. The spatial resolution
%                                           of the height.
%               {'wspatial resolution'}*  - Numeric. The spatial resolution
%                                           of the width.
%               {'sample spacing'}*       - Matrix. The [width, height]
%                                           sample spacing resolution.
%               {'distance per sample'}*  - Matrix. The [row, col] distance
%                                           per spatial sample.
%               {'distance per degree'}*  - Numeric. The distance per
%                                           degree of visual angle.
%               {'degrees per distance'}* - Numeric. *Needs to be added*.
%               {'spatial sampling positions'}*
%                                         - Matrix. Spatial point locations
%               {'hangular resolution'}   - Numeric. The angular degree per
%                                           pixel in height.
%               {'wangular resolution'}   - Numeric. The angular degree per
%                                           pixel in width.
%               {'angular resolution'}    - Matrix. The [height, width]
%                                           angular resolution.
%               {'frequency Support'}*    - Numeric. The frequency
%                                           resolution in cyc/deg or
%                                           lp/Unit, i.e., cycles/{meters,
%                                           mm, microns} ex. oiGet(oi, ...
%                                           'frequencyResolution', 'mm')
%               {'max frequency resolution'}*
%                                         - Numeric. The highest frequency
%                                           oiGet(oi, ...
%                                           'maxFrequencyResolution', 'um')
%               {'frequency support col', 'fsupportx'}*
%                                         - Numeric. The frequency support
%                                           for cols.
%               {'frequency support row', 'fsupporty'}*
%                                         - Numeric. The frequency support
%                                           for rows.
%          Depth
%               {'depthMap'} - Matrix. Pixel wise depth map in meters.
%          Optics information
%               {'optics'}          - Struct. See opticsSet/Get
%               {'optics model'}    - String. Options are: diffraction
%                                     limited, shift invariant, ray trace.
%               {'diffuser method'} - String. Options are: 'skip', 'blur'
%                                     (gaussian), 'birefringent'
%               {'diffuser blur'}   - Numeric. The std of Gaussian blur.
%               {'psfstruct'}       - DEPRECATED. Struct. The entire
%                                     shift-variant PSF structure.
%                  {'sampled rt psf'}    - DEPRECATED. Matrix. Precomputed
%                                          shift-variant psfs.
%                  {'psf sample angles'} - DEPRECATED. Array. Vector of
%                                          sample angle.
%                  {'psf angle step'}    - DEPRECATED. Numeric. The spacing
%                                          between ray trace angle samples.
%                  {'psf image heights'} - DEPRECATED. Array. Sampled image
%                                          heights vector (use optics).
%                  {'raytrace optics name'}
%                                        - DEPRECATED. String. Optics used
%                                          to derive the shift-variant psf.
%                  {'rt psf size'}       - DEPRECATED. Matrix. The [row,
%                                          col] dimensions of the psf.
%          Misc
%               {'rgb image'} - Matrix. The RGB rendering of OI data.
%
% Outputs:
%    val      - VARIES. The value of the requested parameter. See param
%               descriptions above for indications of type.
%
% Optional key/value pairs:
%    Needs to be added.
%

% History:
%	 xx/xx/03       Copyright ImagEval Consultants, LLC, 2003-2015.
%    12/30/17  dhb  Started to put comments into isetbio standard format.
%    01/22/18  dhb  Example runs in clean workspace.
%    03/07/18  jnm  Formatting
%    03/26/18  dhb  Change arg to determineUnits to (..., varargin{1}).
%                   Needs to be like that to work, otherwise a cell array
%                   rather than a string comes back.
%    04/07/18  dhb  Change call to vcGetSelectedObject to vcGetObject. I
%                   did this once before, but somehow it was lost.
%    06/25/19  JNM  Documentation update


% Examples:
%{
    scene = sceneCreate;
    oi = oiCreate;
    oi = oiCompute(oi, scene);

    oiGet(oi, 'rows')
    oiGet(oi, 'wave')
    oiGet(oi, 'optics')
    oiGet(oi, 'area', 'mm')
    oiGet(oi, 'wres', 'microns')
    oiGet(oi, 'angularresolution')
    oiGet(oi, 'distPerSamp', 'mm')

    % Meshgrid of zero-centered (x, y) values
    oiGet(oi, 'spatial support', 'microns');
    oiGet(oi, 'optics off axis method')
    oiGet(oi, 'lens');   % Lens object
%}

if ~exist('parm', 'var') || isempty(parm)
    error('Param must be defined.');
end
val = [];

% See if this is really an optics call. The routine ieParameterOtype will
% do some fragile magic to decide whether want to do a subcall into optics
% or lens.
[oType, parm] = ieParameterOtype(parm);
switch oType
    case 'optics'
        % If optics, then we either return the optics or an optics
        % parameter. We could add another varargin{} level. Or even all.
        optics = oi.optics;
        if isempty(parm), val = optics;
        elseif isempty(varargin), val = opticsGet(optics, parm);
        else, val = opticsGet(optics, parm, varargin{1});
        end
        return;
    case 'lens'
        lens = oi.optics.lens;
        if isempty(parm), val = lens;
        else, val = lens.get(parm, varargin{:});
        end
        return;
end

% It appears to be an oi, so onward.
parm = ieParamFormat(parm);
switch parm
    case 'type'
        val = oi.type;
    case 'name'
        val = oi.name;
    case 'filename'
        val = oi.filename;
    case 'consistency'
        val = oi.consistency;
    case {'rows', 'row', 'nrows', 'nrow'}
        if checkfields(oi, 'data', 'photons')
            val = size(oi.data.photons, 1);
        else
            % disp('Using current scene rows')
            scene = vcGetObject('scene');
            if isempty(scene)
                disp('No scene and no oi. Using 128 rows.');
                val = 128;
            else
                val = sceneGet(scene, 'rows');
            end
        end
    case {'cols', 'col', 'ncols', 'ncol'}
        if checkfields(oi, 'data', 'photons')
            val = size(oi.data.photons, 2);
        else
            % disp('Using current scene cols')
            scene = vcGetObject('scene');
            if isempty(scene)
                disp('No scene and no oi. Using 128 cols.');
                val = 128;
            else
                val = sceneGet(scene, 'cols');
            end
        end
    case 'size'
        val = [oiGet(oi, 'rows'), oiGet(oi, 'cols')];
    case {'lens', 'lenspigment'}
        if checkfields(oi, 'optics', 'lens'), val = oi.optics.lens; end
    case {'samplespacing'}
        % Sample spacing, both height and width
        %    oiGet(oi, 'sample spacing', 'mm')
        % If no OI is yet computed, we pad the scene size as we would have
        % in the oiCompute calculation.
        sz = oiGet(oi, 'size');
        if isempty(sz)
            % Given the current scene, this is what the computed OI size
            % will be.
            scene = vcGetObject('scene');
            sz = sceneGet(scene, 'size');
            if isempty(sz)
                error('No scene or OI');
            else
                padSize = round(sz / 8);
                sz = size(...
                    padarray(zeros(sz(1), sz(2)), padSize, 0, 'both'));
            end
        end
        val = [oiGet(oi, 'width') / sz(2) , oiGet(oi, 'height') / sz(1)];
        if ~isempty(varargin)
            val = val * ieUnitScaleFactor(varargin{1});
        end
    case {'samplesize'}
        % This is the spacing between samples. We expect row and column
        % spacing are equal.
        %    oiGet(oi, 'sample size', 'mm')
        %
        % Not protected from missing oi data as in 'sample spacing'.
        % Should be integrated with that one.
        w = oiGet(oi, 'width');  % Image width in meters
        c = oiGet(oi, 'cols');   % Number of sample columns
        val = w / c;             % M/sample
        if ~isempty(varargin)
            val = val * ieUnitScaleFactor(varargin{1});
        end
    case {'distance', 'imagedistance', 'focalplanedistance'}
        % We are finding the distance from the lens to the focal plane. The
        % focal plane depends on the scene distance. If there is no scene,
        % we assume the user meant a scene at infinity.
        if ~isempty(varargin)
            sDist = varargin{1};
        elseif isfield(oi, 'distance') && oi.distance > 0
            sDist = oi.distance;
        else
            %scene = vcGetSelectedObject('Scene');
            scene = vcGetObject('scene');
            if isempty(scene), sDist = 1e10;
            else, sDist = sceneGet(scene, 'distance');
            end
        end
        val = opticsGet(oiGet(oi, 'optics'), 'imagedistance', sDist);
        if ~isempty(varargin)
            val = val * ieUnitScaleFactor(varargin{1});
        end
    case {'wangular', 'widthangular', 'fov', ...
            'hfov', 'horizontalfieldofview'}
        % oiCompute(oi, scene) assigns the angular field of view to the oi.
        % Normally, an OI would not have an angular field of view; it just
        % reflects the angle of the scene it represents.
        if checkfields(oi, 'wAngular')
            val = oi.wAngular;
        else
            % We use the scene or a default.
            scene = vcGetObject('scene');
            if isempty(scene) % Default scene visual angle.
                disp('Arbitrary oi angle: 10 deg')
                val = 10;
            else  % Use current scene angular width
                val = sceneGet(scene, 'wangular');
            end
        end
    case {'hangular', 'heightangular', 'vfov', 'verticalfieldofview'}
        % We only store the width FOV. We insist that the pixels are square
        h = oiGet(oi, 'height');         % Height in meters
        d = oiGet(oi, 'distance');       % Distance to lens
        val = 2 * atand((0.5 * h) / d);  % Vertical field of view
    case {'dangular', 'diagonalangular', 'diagonalfieldofview'}
        val = sqrt(oiGet(oi, 'wAngular') ^ 2 + oiGet(oi, 'hAngular') ^ 2);
    case 'aspectratio'
        r = oiGet(oi, 'rows');
        c = oiGet(oi, 'cols');
        if isempty(c) || c == 0, disp('No OI'); return;
        else, val = r / c;
        end
    case 'optics'
        % Terms related to the optics
        % This is the large optics structure
        if checkfields(oi, 'optics'), val = oi.optics; end
    case 'opticsmodel'
        if checkfields(oi, 'optics', 'model'), val = oi.optics.model; end

        %{
        % I think this's only used for the shift-variant (raytrace) model.
        % It's not used any more here. It's still retained in ISETCAM.

        % Sometimes we precompute the psf from the optics and store it
        % here. The angle spacing of the precomputation is specified here
    case {'psfstruct'}
        % Entire svPSF structure
        if checkfields(oi, 'psf'), val = oi.psf; end
    case {'svpsf', 'shiftvariantpsf'}
        % Precomputed shift-variant psfs
        if checkfields(oi, 'psf', 'psf'), val = oi.psf.psf; end
    case {'rtpsfsize'}
        % Number of spatial samples in each psf
        if checkfields(oi, 'psf', 'psf')
            val = size(oi.psf.psf{1, 1, 1});
        end
    case {'psfsampleangles'}
        % Vector of sample angle
        if checkfields(oi, 'psf', 'sampAngles')
            val = oi.psf.sampAngles;
        end
    case {'psfanglestep'}
        % Spacing between angles
        if checkfields(oi, 'psf', 'sampAngles')
            val = oi.psf.sampAngles(2) - oi.psf.sampAngles(1);
        end
    case {'psfimageheights'}
        % Vector of sampled image heights
        if checkfields(oi, 'psf', 'imgHeight'), val = oi.psf.imgHeight; end
        if ~isempty(varargin)
            val = val * ieUnitScaleFactor(varargin{1});
        end
    case {'psfopticsname', 'raytraceopticsname'}
        % Optics data are derived from
        if checkfields(oi, 'psf', 'opticsName')
            val = oi.psf.opticsName;
        end
    case 'psfwavelength'
        % Wavelengths for this calculation. Should match the optics, I
        % think.  Not sure why it is duplicated.
        if checkfields(oi, 'psf', 'wavelength')
            val = oi.psf.wavelength;
        end
        %}
        % optical diffuser properties
  case {'diffusermethod'}
      % 0 - skip, 1 - gauss blur, 2 - birefringent
      if checkfields(oi, 'diffuser', 'method')
          val = oi.diffuser.method;
      end
    case {'diffuserblur'}
      if checkfields(oi, 'diffuser', 'blur'), val = oi.diffuser.blur; end
      if ~isempty(varargin)
          val = val * ieUnitScaleFactor(varargin{1});
      end
    case 'data'
       if checkfields(oi, 'data'), val = oi.data; end
   case {'photons'}
       % Read photon data.
       % Data are returned as doubles.
       if checkfields(oi, 'data', 'photons')
           if isempty(varargin)
               % allPhotons = oiGet(oi, 'photons')
               val = double(oi.data.photons);
           else
               % waveBandPhotons = oiGet(oi, 'photons', 500)
               idx = ieFindWaveIndex(oiGet(oi, 'wave'), varargin{1});
               val = double(oi.data.photons(:, :, idx));
           end
       end
    case 'roiphotons'
       if isempty(varargin), error('ROI required')
       else, roiLocs = varargin{1};
       end
       val = vcGetROIData(oi, roiLocs, 'photons');
   case 'roimeanphotons'
       if isempty(varargin), error('ROI required')
       else, roiLocs = varargin{1};
       end
       val = oiGet(oi, 'roiphotons', roiLocs);
       val = mean(val, 1);
    case {'photonsnoise', 'photonswithnoise'}
        % pn = oiGet(oi, 'photons noise');
        % The current photons are the mean.
        % This returns the mean photons plus Poisson noise
        %
        % To compute noise, we need a photon count, not rate. We choose the
        % pixel area and a 50 msec integration time, somewhat arbitrarily.
        sz = oiGet(oi, 'sample spacing');
        t = 0.050;
        photons = oiGet(oi, 'photons');
        val = oiPhotonNoise(photons * sz(1) * sz(2) * t);
    case {'energynoise', 'energywithnoise'}
        % en = oiGet(oi, 'energy noise');
        % Return mean energy plus Poisson photon noise
        sz = oiGet(oi, 'sample spacing');
        t = 0.050;
        photons = oiGet(oi, 'photons');
        val = oiPhotonNoise(photons * sz(1) * sz(2) * t);
        wave = oiGet(oi, 'wave');
        val = Quanta2Energy(wave(:), val);
    case {'datamax', 'dmax'}
        % return data max, not for compression anymore
        if checkfields(oi, 'data', 'photons')
            val = max(oi.data.photons(:));
        end
    case {'datamin', 'dmin'}
        % return data min, not for compression anymore
        if checkfields(oi, 'data', 'photons')
            val = min(oi.data.photons(:));
        end
    case {'bitdepth', 'compressbitdepth'}
        if checkfields(oi, 'data', 'bitDepth'), val = oi.data.bitDepth; end
    case 'energy'
        % Compute energy from photons
        if checkfields(oi, 'data', 'photons')
            wave = oiGet(oi, 'wave');
            val = Quanta2Energy(wave, oi.data.photons);
        end
    case 'roienergy'
       if isempty(varargin), error('ROI required')
       else, roiLocs = varargin{1};
       end
       val = vcGetROIData(oi, roiLocs, 'energy');
    case 'roimeanenergy'
       if isempty(varargin), error('ROI required')
       else, roiLocs = varargin{1};
       end
       val = oiGet(oi, 'roienergy', roiLocs);
       val = mean(val, 1);
    case {'meanilluminance', 'meanillum'}
        % Get / compute mean illuminance
        % Always update the mean. Cheap to do. Forces synchronization.
        if notDefined('oi.data.illuminance') || ...
                isempty(oi.data.illuminance)
            oi.data.illuminance = oiCalculateIlluminance(oi);
        end
        val = mean(oi.data.illuminance(:));
    case {'illuminance', 'illum'}
        if isempty(oi.data.illuminance)
            % calculate and store
            oi.data.illuminance = oiCalculateIlluminance(oi);
        end
        val = oi.data.illuminance;
    case {'xyz', 'dataxyz'}
        % oiGet(oi, 'xyz');
        % RGB array of oi XYZ values. These are returned as an RGB format
        % at the spatial sampling grid of the optical image.
        photons = oiGet(oi, 'photons');
        wave = oiGet(oi, 'wave');
        val = ieXYZFromEnergy(Quanta2Energy(wave, photons), wave);
        % sz = oiGet(oi, 'size');
        % val = XW2RGBFormat(val, sz(1), sz(2));
    case {'lms', 'datalms', 'cone'}
        % oiGet(oi, 'lms');
        % RGB (3D array) of oi Stockman LMS values
        energy = oiGet(oi, 'energy');
        wave = oiGet(oi, 'wave');
        S = ieReadSpectra('stockman', wave);
        if numel(wave) > 1
            dWave = wave(2) - wave(1);
        else
            dWave = 10;
            disp('10 nm bandwidth assumed');
        end
        [energy, r, c] = RGB2XWFormat(energy);
        val = energy * S * dWave;
        val = XW2RGBFormat(val, r, c);
    case {'spectrum', 'wavespectrum'}
        if isfield(oi, 'spectrum'), val = oi.spectrum; end
    case 'binwidth'
        wave = oiGet(oi, 'wave');
        if length(wave) > 1, val = wave(2) - wave(1); else, val = 1; end
    case {'datawave', 'photonswave', 'photonswavelength', ...
            'wave', 'wavelength'}
        % oiGet(oi, 'wave')
        % There are a number of different
        % differ from the optics spectrum. There should only be one, and
        % it should probably be part of the optics. To smooth the
        % transition to that wonderful day, we return the optics spectrum
        % if there is no oi spectrum. Always a column vector, even if
        % people stick it in the wrong way.
        if checkfields(oi, 'spectrum', 'wave')
            val = oi.spectrum.wave(:);
        % elseif checkfields(oi, 'optics', 'spectrum', 'wave'),
        %     val = oi.optics.spectrum.wave(:);
        % elseif checkfields(oi, 'optics', 'OTF', 'wave')
        %     val = oi.optics.OTF.wave(:);
        % elseif checkfields(oi, 'optics', 'rayTrace', 'psf', 'wavelength')
        %     val = oi.optics.rayTrace.psf.wavelength(:);
        end
    case {'nwave', 'nwaves'}
        % Refers to data wave.
        % oiGet(oi, 'n wave')
        val = length(oiGet(oi, 'data wave'));
    case 'height'
        % Height in meters is default
        % oiGet(oi, 'height', 'microns')
        val = oiGet(oi, 'sampleSize') * oiGet(oi, 'rows');
        if ~isempty(varargin)
            val = val * ieUnitScaleFactor(varargin{1});
        end
    case {'width'}
        % Width in meters is default - We need to handle 'skip' case
        d = oiGet(oi, 'focalPlaneDistance');  % Distance from lens to image
        fov = oiGet(oi, 'wangular');          % FOV (horizontal, width)
        % ieRad2deg(2 * atan((0.5 * width) / imageDistance)) = fov
        val = 2 * d * tand(fov / 2);
        if ~isempty(varargin)
            val = val * ieUnitScaleFactor(varargin{1});
        end
    case {'diagonal', 'diagonalsize'}
        val = sqrt(oiGet(oi, 'height') ^ 2 + oiGet(oi, 'width') ^ 2);
        if ~isempty(varargin)
            val = val * ieUnitScaleFactor(varargin{1});
        end
    case {'heightwidth', 'heightandwidth'}
        val(1) = oiGet(oi, 'height');
        val(2) = oiGet(oi, 'width');
        if ~isempty(varargin)
            val = val * ieUnitScaleFactor(varargin{1});
        end
    case {'area', 'areameterssquared'}
        % oiGet(oi, 'area')        % square meters
        % oiGet(oi, 'area', 'mm')  % square millimeters
        val = oiGet(oi, 'height') * oiGet(oi, 'width');
        if ~isempty(varargin)
            val = val * ieUnitScaleFactor(varargin{1}) ^ 2;
        end
    case {'centerpixel', 'centerpoint'}
        val = [oiGet(oi, 'rows'), oiGet(oi, 'cols')];
        val = floor(val / 2) + 1;
    case {'hspatialresolution', 'heightspatialresolution', 'hres'}
        % Resolution parameters
        % Size in distance per pixel, default is meters per pixel
        % oiGet(oi, 'hres', 'microns') is acceptable syntax.
        h = oiGet(oi, 'height');
        r = oiGet(oi, 'rows');
        if isempty(r) && strcmp(oiGet(oi, 'type'), 'opticalimage')
            % For optical images we return a default based on the scene.
            % This is used when no optical image has been calculated.
            scene = vcGetObject('scene');
            if isempty(scene)
                disp('No scene or oi. Using 128 rows');
                r = 128; % Make something up
            else
                r = oiGet(scene, 'rows');
            end
        end
        val = h / r;
        if ~isempty(varargin)
            val = val * ieUnitScaleFactor(varargin{1});
        end
    case {'wspatialresolution', 'widthspatialresolution', 'wres'}
        % Size in m per pixel is default.
        % oiGet(oi, 'wres', 'microns')
        % oiGet(oi, 'wres')
        w = oiGet(oi, 'width');
        c = oiGet(oi, 'cols');
        if isempty(c)
            % For optical images we return a default based on the scene.
            % This is used when no optical image has been calculated.
            scene = vcGetObject('scene');
            if isempty(scene)
                disp('No scene or oi. Using 128 cols');
                c = 128;
            else
                c = oiGet(scene, 'cols');
            end
        end
        val = w / c;
        if ~isempty(varargin)
            val = val * ieUnitScaleFactor(varargin{1});
        end
    case {'spatialresolution', 'distancepersample', 'distpersamp'}
        % oiGet(oi, 'distPerSamp', 'mm')
        val = [oiGet(oi, 'hspatialresolution'), ...
               oiGet(oi, 'wspatialresolution')];
        if ~isempty(varargin)
            val = val * ieUnitScaleFactor(varargin{1});
        end
    case {'distperdeg', 'distanceperdegree'}
        % This routine should call
        % opticsGet(optics, 'dist per deg', unit) rather than compute it
        % here.
        val = oiGet(oi, 'width') / oiGet(oi, 'fov');
        if ~isempty(varargin)
            val = val * ieUnitScaleFactor(varargin{1});
        end
    case {'degreesperdistance', 'degperdist'}
        % oiGet(oi, 'degrees per distance')
        % oiGet(oi, 'degPerDist', 'micron')
        %
        % We should probably call:
        %   opticsGet(optics, 'dist per deg', unit)
        % which is preferable to this call.
        % if isempty(varargin), units = 'm'; else, units = varargin{1}; end
        units = determineUnits('m', varargin{1});
        val = oiGet(oi, 'distance per degree', units);  % meters
        val = 1 / val;
        % val = oiGet(oi, 'fov')/oiGet(oi, 'width');
        % if ~isempty(varargin)
        %    val = val / ieUnitScaleFactor(varargin{1});
        % end
    case {'spatialsupport', 'spatialsamplingpositions'}
        % oiGet(oi, 'spatialsupport', 'microns')
        % Spatial locations of points in meters
        sSupport = oiSpatialSupport(oi);
        [xSupport, ySupport] = meshgrid(sSupport.x, sSupport.y);
        val(:, :, 1) = xSupport;
        val(:, :, 2) = ySupport;
        if ~isempty(varargin)
            val = val * ieUnitScaleFactor(varargin{1});
        end
    case {'angularsupport', 'angularsamplingpositions'}
        % Angular values of sample points
        % Units can be deg (default), min or sec and should include radians
        degPerSamp = oiGet(oi, 'angular resolution');
        sz = oiGet(oi, 'size');
        tmp = degPerSamp(2) * (1:sz(2));
        aSupport.x = tmp - mean(tmp(:));
        tmp = degPerSamp(1) * (1:sz(1));
        aSupport.y = tmp - mean(tmp(:));
        % Deal with units
        if ~isempty(varargin)
            unit = lower(varargin{1});
            switch unit
                case {'deg'}
                    % Default
                case 'min'
                    aSupport.x = aSupport.x * 60;
                    aSupport.y = aSupport.y * 60;
                case 'sec'
                    aSupport.x = aSupport.x * 60 * 60;
                    aSupport.y = aSupport.y * 60 * 60;
                case 'radians'
                    aSupport.x = ieDeg2rad(aSupport.x);
                    aSupport.y = ieDeg2rad(aSupport.y);
                otherwise
                    error('Unknown angular unit %s\n', unit);
            end
        end
        [xSupport, ySupport] = meshgrid(aSupport.x, aSupport.y);
        val(:, :, 1) = xSupport;
        val(:, :, 2) = ySupport;
    case {'hangularresolution', 'heightangularresolution'}
        % Angular degree per pixel -- in degrees
        val = 2 * atand((oiGet(oi, 'hspatialResolution') / ...
            oiGet(oi, 'distance')) / 2);
    case {'wangularresolution', 'widthangularresolution'}
        % Angle (degree) per pixel -- width
        val = 2 * atand((oiGet(oi, 'wspatialResolution') / ...
            oiGet(oi, 'distance')) / 2);
    case {'angularresolution', 'degperpixel', 'degpersample', ...
            'degreepersample', 'degreeperpixel'}
        % Height and width
        val = [oiGet(oi, 'hangularresolution'), ...
            oiGet(oi, 'wangularresolution')];
    case {'frequencyresolution', 'freqres'}
        % Default is cycles per degree
        % val = oiGet(oi, 'frequencyResolution', units);

        % if isempty(varargin), units = 'cyclesPerDegree';
        % else, units = varargin{1};
        % end
        units = determineUnits('cyclesPerDegree', varargin{1});
        val = oiFrequencySupport(oi, units);
    case {'maxfrequencyresolution', 'maxfreqres'}
        % Default is cycles/deg. By using
        % oiGet(oi, 'maxfreqres', units)
        % you can get cycles/{meters, mm, microns}

        % if isempty(varargin), units = 'cyclesPerDegree';
        % else, units = varargin{1};
        % end
        units = determineUnits('cyclesPerDegree', varargin{1});
        % val = oiFrequencySupport(oi, units);
        if isempty(varargin), units = []; end
        fR = oiGet(oi, 'frequencyResolution', units);
        val = max(max(fR.fx), max(fR.fy));
    case {'frequencysupport', 'fsupportxy', 'fsupport2d', 'fsupport'}
        % val = oiGet(oi, 'frequency support', units);

        % if isempty(varargin), units = 'cyclesPerDegree';
        % else, units = varargin{1};
        % end
        units = determineUnits('cyclesPerDegree', varargin{1});
        fResolution = oiGet(oi, 'frequencyresolution', units);
        [xSupport, ySupport] = meshgrid(fResolution.fx, fResolution.fy);
        val(:, :, 1) = xSupport;
        val(:, :, 2) = ySupport;
    case {'frequencysupportcol', 'fsupportx'}
        % val = oiGet(oi, 'frequency support col', units);

        % if isempty(varargin), units = 'cyclesPerDegree';
        % else, units = varargin{1};
        % end
        units = determineUnits('cyclesPerDegree', varargin{1});
        fResolution = oiGet(oi, 'frequencyresolution', units);
        l = find(abs(fResolution.fx) == 0);
        val = fResolution.fx(l:end);
    case {'frequencysupportrow', 'fsupporty'}
        % val = oiGet(oi, 'frequency support row', units);

        % if isempty(varargin), units = 'cyclesPerDegree';
        % else, units = varargin{1};
        % end
        units = determineUnits('cyclesPerDegree', varargin{1});
        fResolution = oiGet(oi, 'frequencyresolution', units);
        l = find(abs(fResolution.fy) == 0);
        val = fResolution.fy(l:end);
    case {'customcomputemethod', 'oicompute', 'oicomputemethod', ...
            'oimethod'}
        % Computational methods -- About to be obsolete and managed by the
        % optics model information in the optics structure.
        if checkfields(oi, 'customMethod'), val = oi.customMethod; end
    case {'customcompute', 'booleancustomcompute'}
        % 1 or 0
        if checkfields(oi, 'customCompute'), val = oi.customCompute;
        else, val = 0;
        end
    case {'rgb', 'rgbimage'}
        % Visual information
        %
        % Get the rgb image shown in the oi window. If the displayFlag is
        % set, show the image too.
        %     rgb = oiGet(oi, 'rgb image', gamma = 1, displayFlag = -1);
        if checkfields(oi, 'data', 'photons')
            % Get the photons to later render the rgb image
            % Save a lot of memory by direct access.
            photons = oi.data.photons;
        else
            % No photons, so return val as empty
            return;
        end
        % Deal with display gamma
        if isempty(varargin)
            oiW = ieSessionGet('oi window handle');
            if isempty(oiW), gam = 1;
            else, gam = str2double(get(oiW.editGamma, 'string'));
            end
        else
            gam = varargin{1};
        end
        wList = oiGet(oi, 'wave');
        [row, col, ~] = size(photons);
        displayFlag = -1;  % Default is do not display, just compute
        if length(varargin) > 1, displayFlag = varargin{2}; end
        % Compute and show (or not)
        val = imageSPD(photons, wList, gam, row, col, displayFlag);
    case {'depthmap'}
        % The depth map in the OI domain indicates locations where there
        % were legitimate scene data for computing the OI data. Other
        % regions are 'extrapolated' by sceneDepthRange to keep the
        % calculations correct. But they don't always correspond to the
        % original data. When there is no depthMap in the scene, these are
        % all logically '1' (true).
        if checkfields(oi, 'depthMap'), val = oi.depthMap; end
    case {'logicaldepthmap'}
        % Boolean values indicating locations computed from
        val = logical(oiGet(oi, 'depth map'));
    case {'depthmapphotons'}
        % Get the photons that are within the specified depth map region
        error('Not yet implemented. See oiCombineDepth for algorithm');
    otherwise
        error(['Unknown parameter: ', parm]);
 end

end