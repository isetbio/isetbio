function val = displayGet(d, parm, varargin)
% Get display parameters and derived properties
%
% Syntax:
%   val = displayGet(d, parm, varargin)
%
% Description:
%    Get the display parameters and derived properties for a display.
%
% Inputs:
%    d        - Struct. The display structure
%    parm     - String. A string describing the parameter name. Since
%               there are so many options, they are separated out into
%               sections below. Will also contain return type when able.
%        Basic parameters
%            {'type'}                 - String. Always 'display'
%            {'name'}                 - String. Which specific display
%            {'is emissive'}          - Boolean. true (emissive) or
%                                       false (reflective)
%            {'main image'}           - Matrix. image for the main
%                                       display window
%
%        Transduction
%            {'gamma table'}          - Matrix. nLevels x nPrimaries
%            {'inverse gamma', [nSteps]}
%                                     - invert gamma table
%                                       See ieLUTInvert
%            {'dacsize'}              - Numeric. The number of bits
%                                       (log2(nSamples))
%            {'nlevels'}              - Numeric. number of levels
%            {'levels'}               - List. list of levels
%
%        SPD calculations
%            {'wave'}                 - Vector. wavelength samples in nm
%            {'nwave'}                - Numeric. number of wave samples
%            {'spd primaries'}        - Matrix. nWave x nPrimaries
%                                       matrix, in energy units
%            {'white spd'}            - SPD. The white point spectral
%                                       power distribution
%            {'black spd'}            - SPD. Spd when display is black
%            {'nprimaries'}           - Numeric. number of primaries
%
%        Color conversion and metric
%            {'rgb2xyz'}              - Linear rgb to CIE (1931) XYZ
%            {'rgb2lms'}              - Linear rgb to Stockman cones
%                                       (row format)
%            {'lms2rgb'}              - Stockman cones to linear RGB
%                                       (for cone isolation)
%            {'white xyz'}            - The Display white point in
%                                       CIE (1931) XYZ
%            {'black xyz'}            - This includes the stray light
%                                       from backlight
%            {'white xy'}
%            {'white lms'}
%            {'primaries xyz'}
%            {'primaries xy'}
%            {'peak luminance'}
%            {'dark luminance'}
%            {'peak contrast'}
%
%        Spatial parameters
%            {'dpi', 'ppi'}           - Numeric. dots per inch
%            {'size'}                 - Vector. 1x2 vector of display
%                                       size in meters [h, v]
%            {'meters per dot'}
%            {'dots per meter'}
%            {'dots per deg'}         - Numeric. The number of dots per
%                                       degree visual angle
%            {'viewing distance'}     - Numeric. VD in meters
%
%        Subpixel structure
%            {'dixel'}                - Struct. The dixel structure
%                                       describing repeating unit
%            {'pixels per dixel'}     - Numeric. number of pixels in one
%                                       repeating unit
%            {'dixel size'}           - Numeric. The number of samples
%                                       in one dixel
%            {'dixel image'}          - Matrix. image in dixel panel
%            {'dixel control map'}    - The control map, describing
%                                       which of the regions are
%                                       individually addressable
%            {'peak spd'}             - peak spd for each primary
%            {'oversample'}           - The up-scale factor for the
%                                       subpixel rendering
%            {'sample spacing'}       - Numeric. spacing between samples
%            {'fill factor'}          - The fill factor of each primary
%            {'render function'}      - function handle used to convert
%                                       rgb input to the corresponding
%                                       dixel image
%    varargin - (Optional) Additional argument(s) as necessary
%
% Outputs:
%    val      - VARIES. The value of the desired parameter.
%
% Optional key/value pairs:
%    **NEEDS TO BE FILLED OUT**
%

% History:
%    xx/xx/14  HJ/BW  ISETBIO TEAM, Copyright 2014
%    06/25/15  dhb    Added get of display size
%    05/16/18  jnm    Formatting

% Examples:
%{
    d = displayCreate;
    w = displayGet(d, 'wave');
    p = displayGet(d, 'spd');
    vcNewGraphWin;
    plot(w, p);
    set(gca, 'ylim', [-.1 1.1])

    chromaticityPlot(displayGet(d, 'white xy'))

    % vci = vcimageCreate('test', [], d);
    % plotDisplayGamut(vci)
%}

%% Check parameters
if notDefined('parm'), error('Parameter not found.');  end

% Default is empty when the parameter is not yet defined.
val = [];

parm = ieParamFormat(parm);

%% Do the analysis
switch parm
    case {'name'}
        val = d.name;

    case {'type'}
        % Type should always be 'display'
        val = d.type;

    case {'gtable', 'dv2intensity', 'gamma', 'gammatable'}
        if isfield(d, 'gamma'), val = d.gamma; end

    case {'inversegamma', 'inversegammatable'}
        if isfield(d, 'gamma')
            % Optional nSteps arg for inverse gamma table
            if (isempty(varargin))
                val = ieLUTInvert(d.gamma);
            else
                val = ieLUTInvert(d.gamma, varargin{1});
            end
        end

    case {'isemissive'}
        val = true;
        if isfield(d, 'isEmissive'), val = d.isEmissive; end

    case {'mainimage', 'mainimg'}
        % Main image in display window
        val = d.mainimage;

    case {'bits', 'dacsize'}
        % color bit depths, e.g. 8 bit / 10 bit
        % This is computed from size of gamma table
        gTable = displayGet(d, 'gTable');
        assert(ismatrix(gTable), 'Bit depth of display unkown');
        val = round(log2(size(gTable, 1)));

    case {'nlevels'}
        % Number of levels
        val = 2 ^ displayGet(d, 'bits');

    case {'levels'}
        % List of the levels, e.g. 0~255
        val = 1:displayGet(d, 'nlevels') - 1;

    case {'wave', 'wavelength'}  % SPD calculations in nanometers
        % For compatibility with PTB. We might change .wave to
        % .wavelengths.
        if checkfields(d, 'wave'), val = d.wave(:);
        elseif checkfields(d, 'wavelengths'), val = d.wavelengths(:);
        end

    case {'binwidth'}
        wave = displayGet(d, 'wave');
        if length(wave) > 1, val = wave(2) - wave(1); end

    case {'nwave'}
        val = length(displayGet(d, 'wave'));

    case {'nprimaries'}
        % SPD is always nWave by nPrimaries
        spd = displayGet(d, 'spd');
        val = size(spd, 2);

    case {'spd', 'spdprimaries'}
        % Units are energy (watts/....)
        % displayGet(dsp, 'spd');
        % displayGet(d, 'spd', wave);
        %
        % Always make sure the spd has rows equal to number of
        % wavelength samples. The PTB uses spectra rather than spd. This
        % hack makes it compatible. Or, we could convert displayCreate
        % from spd to spectra some day.
        if checkfields(d, 'spd'), val = d.spd;
        elseif checkfields(d, 'spectra'), val = d.spectra;
        end

        % Sometimes users put the data in transposed, sigh. I am one of
        % those users.
        nWave = displayGet(d, 'nwave');
        if size(val, 1) ~= nWave, val = val'; end

        % Interpolate for alternate wavelength, if requested
        if ~isempty(varargin)
            % Wave interpolation
            wavelength = displayGet(d, 'wave');
            wave = varargin{1};
            val = interp1(wavelength(:), val, wave(:), 'linear', 0);
        end

    case {'whitespd'}
        % SPD when all the primaries are at peak, this is the energy
        if ~isempty(varargin), wave = varargin{1};
        else, wave = displayGet(d, 'wave');
        end
        e = displayGet(d, 'spd', wave);
        val = sum(e, 2);

    case {'rgb2xyz'}  % Color conversion
        % rgb2xyz = displayGet(dsp, 'rgb2xyz', wave)
        % RGB as a column vector mapped to XYZ column
        %  x(:)' = r(:)' * rgb2xyz
        % Hence, imageLinearTransform(img, rgb2xyz)
        % should work
        wave = displayGet(d, 'wave');
        spd = displayGet(d, 'spd', wave);  % spd in energy
        val = ieXYZFromEnergy(spd', wave);

    case {'rgb2lms'}
        % rgb2lms = displayGet(dsp, 'rgb2lms')
        % rgb2lms = displayGet(dsp, 'rgb2lms', wave)
        %
        % RGB as a row vector is mapped to an LMS row
        %
        %     c(:)' = r(:)' * rgb2lms
        %
        % This matrix format is used here
        %
        %   imageLinearTransform(img, rgb2lms)
        %
        wave = displayGet(d, 'wave');

        % The Stockman Energy Fundamentals include the photopigment,
        % default lens and default macular pigment (BW, I think).
        coneFile = fullfile(isetbioDataPath, 'human', 'stockman');
        cones = ieReadSpectra(coneFile, wave);  % plot(wave, spCones)
        spd = displayGet(d, 'spd', wave);       % plot(wave, displaySPD)
        val = cones' * spd;
        val = val';

    case {'lms2rgb'}
        % Linear rgb to Stockman cone coordinates
        % [R, G, B] = [L, M, S] * lm2rgb
        %
        val = inv(displayGet(d, 'rgb2lms'));

    case {'whitexyz', 'whitepoint'}
        % displayGet(dsp, 'white xyz', wave)
        e = displayGet(d, 'white spd');
        if isempty(varargin), wave = displayGet(d, 'wave');
        else, wave = varargin{1};
        end
        % Energy needs to be XW format, so a row vector
        val = ieXYZFromEnergy(e', wave);

    case {'peakluminance'}
        % Luminance of the white point in cd/m2
        % displayGet(dsp, 'peak luminance')
        whiteXYZ = displayGet(d, 'white xyz');
        val = whiteXYZ(2);

    case {'whitexy'}
        val = chromaticity(displayGet(d, 'white xyz'));

    case {'primariesxyz'}
        spd = displayGet(d, 'spd primaries');
        wave = displayGet(d, 'wave');
        val = ieXYZFromEnergy(spd', wave);

    case {'primariesxy'}
        xyz = displayGet(d, 'primaries xyz');
        val = chromaticity(xyz);

    case {'whitelms'}
        % displayGet(dsp, 'white lms')
        rgb2lms = displayGet(d, 'rgb2lms');
        % Sent back in XW format, so a row vector
        val = sum(rgb2lms);

    case {'dpi', 'ppi'}  % Spatial parameters
        if checkfields(d, 'dpi'), val = d.dpi; else, val = 96; end

    case {'size'}
        if checkfields(d, 'size'), val = d.size;
        else, val = [1024 / 768 * 0.3, 0.3];
        end

    case {'metersperdot'}
        % displayGet(dsp, 'meters per dot', 'm')
        % displayGet(dsp, 'meters per dot', 'mm')
        % Useful for calculating image size in meters
        dpi = displayGet(d, 'dpi');
        ipm = 1 / .0254;   % Inch per meter
        dpm = dpi * ipm;   % Dots per meter
        val = 1 / dpm;     % meters per dot
        if ~isempty(varargin)
            val = val * ieUnitScaleFactor(varargin{1});
        end

    case {'dotspermeter'}
        % displayGet(dsp, 'dots per meter', 'm')
        mpd = displayGet(d, 'meters per dot');
        val = 1 / mpd;
        if ~isempty(varargin)
            val = val * ieUnitScaleFactor(varargin{1});
        end

    case {'dotsperdeg', 'sampperdeg'}
        % Samples per deg
        % displayGet(d, 'dots per deg')
        mpd = displayGet(d, 'meters per dot');
        dist = displayGet(d, 'Viewing Distance');  % Meters
        degPerPixel = atand(mpd / dist);
        val = round(1 / degPerPixel);

    case {'degperpixel', 'degperdot'}
        % degrees per pixel
        % displayGet(d, 'deg per dot')
        mpd = displayGet(d, 'meters per dot');
        dist = displayGet(d, 'Viewing Distance');  % Meters
        val = atand(mpd / dist);

    case {'viewingdistance', 'distance'}
        % Viewing distance in meters
        if checkfields(d, 'dist')
            val = d.dist;
        else
            % Default viewing distance in meters, 19 inches
            val = 0.5;
        end

    case {'refreshrate'}
        % display refresh rate
        if isfield(d, 'refreshRate'), val = d.refreshRate; end

    case {'dixel'}  % Dixel (subpixel) information
        % The whole dixel structure
        % displayGet(d, 'dixel')
        if isfield(d, 'dixel'), val = d.dixel; end

    case {'dixelsize'}
        % number of samples in one dixel
        % displayGet(d, 'dixel size')
        dixel_image = displayGet(d, 'dixel image');
        val = size(dixel_image);
        val = val(1:2);

    case {'oversample', 'osample'}
        % Number of subpixel samples per pixel
        % displayGet(d, 'over sample')
        sz = displayGet(d, 'dixel size');
        val = sz ./ displayGet(d, 'pixels per dixel');

    case {'samplespacing'}
        % spacing between psf samples
        % displayGet(d, 'sample sampling', units)
        val = displayGet(d, 'metersperdot') ...
            ./ displayGet(d, 'dixel size');

        % adjust for the number of pixels in one dixel
        val = val .* displayGet(d, 'pixels per dixel');

        if ~isempty(varargin)
            val = val * ieUnitScaleFactor(varargin{1});
        end

    case {'fillfactor', 'fillingfactor', 'subpixelfilling'}
        % Fill factor of subpixle for each primary
        % displayGet(d, 'fill factor')
        dixel_image = displayGet(d, 'dixel image');
        [r, c, ~] = size(dixel_image);
        dixel_image = dixel_image ...
            ./ repmat(max(max(dixel_image)), [r c]);
        dixel_image = dixel_image > 0.2;
        val = sum(sum(dixel_image)) / r / c;
        val = val(:);

    case {'subpixelspd'}
        % spectral power distribution for subpixels
        %
        % This is the real subpixel spd, not the spatial averaged one
        % To get the spd for the whole pixel, use displayGet(d, 'spd')
        % instead
        spd = displayGet(d, 'spd');
        ff = displayGet(d, 'filling factor');
        val = spd ./ repmat(ff(:)', [size(spd, 1) 1]);

    case {'pixelsperdixel'}
        % number of pixels per dixel
        % returns number of pixels in one block (unit repeated pattern)
        % displayGet(d, 'pixels per dixel')
        %
        % The field indicates how many pixels (defined as independent
        % addressable (R, G, B, etc) tuple) in one repeating pattern. In
        % most cases, this field is [1 1], meaning that one dixel
        % contains one pixel. For some displays (say samsung s-strip
        % design), one repeating pattern could contain four independent
        % addressable pixels and in that case pixelsperdixel is [2 2].
        % When HJ built that structure, he assumed that the area
        % occupied by the R, G, B are the same (#R = #G = #B). This may
        % be violated in some displays. Might change this structure when
        % we want to modify it next time.
        if checkfields(d, 'dixel', 'nPixels')
            val = d.dixel.nPixels;
        else
            dixel_control = displayGet(d, 'dixel control map');
            val = max(dixel_control(:));
        end

    case {'dixelintensitymap', 'dixelimage'}
        % dixel intensity map
        % This field specify the intensity (scale factor) at each sample
        % point in dixel
        %
        % displayGet(d, 'dixel intensity map')
        dixel = displayGet(d, 'dixel');
        if isempty(dixel), error('dixel structure not exist'); end
        if isfield(dixel, 'intensitymap')
            val = dixel.intensitymap;
        end

        % adjust the size of the intensity map if required
        if ~isempty(varargin)
            sz = varargin{1};
            if isscalar(sz), sz = [sz sz]; end

            % resize the intensity map
            val = imresize(val, sz);

            % crop the intensity map and make it non-negative
            val(val < 0) = 0;


        % scale the intensity map
            scale = prod(sz) ./ sum(sum(val));
            val = bsxfun(@times, val, scale);
        end

    case {'dixelcontrolmap'}
        % dixel control map
        % This field specifies which region(s) in one dixel is
        % individually addressable
        %
        % The control map contains integer values from 1 ~ n, its value
        % indicates which control group (actual pixel) it belongs to
        %
        % displayGet(d, 'dixel control map')
        dixel = displayGet(d, 'dixel');
        if isempty(dixel), error('dixel structure not exist'); end
        if isfield(dixel, 'controlmap'), val = dixel.controlmap; end

        % adjust the size of the control map if required
        if ~isempty(varargin)
            sz = varargin{1};
            if isscalar(sz), sz = [sz sz]; end
            val = imresize(val, sz, 'nearest'); % resize
        end

    case {'renderfunction'}
        % render function
        % returns user defined subpixel render function handle. If user
        % does not specify this render function, return empty
        %
        % displayGet(d, 'render function')
        if checkfields(d, 'dixel', 'renderFunc')
            val = d.dixel.renderFunc;
        end

    case {'contrast', 'peakcontrast'}
        % peak contrast
        % returns the black/white contrast of the display
        %
        % displayGet(d, 'peak contrast')
        peakLum = displayGet(d, 'peak luminance');
        darkLum = displayGet(d, 'dark luminance');
        val = peakLum / darkLum;

    case {'blackspd', 'blackradiance', 'ambientspd'}
        % black radiance
        % computes dark spd (radiance) of the display in units of energy
        % (watts / ...)
        %
        % displayGet(d, 'black radiance', [wave])
        wave = displayGet(d, 'wave');
        if isfield(d, 'ambient')
            val = d.ambient;
        else
           %  warning(['black (ambient) SPD is not set for display, ' ...
           %      'return 0']);
            val = zeros(size(wave));
        end

        if ~isempty(varargin)
            newWave = varargin{1};
            val = interp1(wave, val, newWave, 'linear', 0);
        end

    case {'darkluminance', 'blackluminance'}
        % dark luminance
        % returns luminance of display when all pixels are turned off
        %
        % displayGet(d, 'dark luminance')
        blackSpd = displayGet(d, 'black spd');
        blackXYZ = ieXYZFromEnergy(blackSpd', displayGet(d, 'wave'));
        val = blackXYZ(2);
    otherwise
        error('Unknown parameter %s\n', parm);
end

end