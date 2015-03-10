function val = opticsGet(optics,parm,varargin)
% Get optics parameters (some values are stored but many are computed)
%
%      val = opticsGet(optics,parm,varargin)
%
% The optics parameters are stored in two different groups of fields.
%
% There are three types of optics models, of increasing complexity. The
% method in use can be selected by the popup menu in the optics window.
% The method is selected by the parameter opticsSet(optics,'model',parameter);
%
% (1) By default, we use a diffraction limited calculation with the
% f-number and focal length determined from the window interface.  Other
% parameters are derived from these few values.  In this case, the optical
% image is computed using oiCompute. To set this method, use
% opticsSet(optics,'model','diffractionLimited');
%
% (2) We can use shift-invariant calculation based on a numerically-defined
% OTF that is wavelength-dependent (but shift-invariant), stored in the
% optics.OTF structure.  When using this method, the user must supply the
% optics structure containing the OTF and other parameters (focal length,
% aperture, and so forth).  Examples are stored in data\optics
% directory.  In this case, the optical image is computed using the
% function opticsSICompute.  To set this method use
% opticsSet(optics,'model','shiftInvariant');
%
% (3) We allow a shift-variant optics model based on ray tracing
% calculations.  In this case, the optics model includes geometric
% distortions, relative illumination and spatially-varying blurring
% calculated using ray tracing programs, such as Code V and Zemax.  The
% geometric distortions, relative illumination, and point spread parameters
% are stored in optics.rayTrace structure. In this case, the optical image
% is computed using the function opticsRayTrace.  To set this method use
% opticsSet(optics,'model','rayTrace')
%
% About the calculations and PSF
%  We store the OTF data with DC at (1,1).  This is true throughout ISET.
%  To understand the implications for certain calculations see the script
%  and tutorial in s_FFTinMatlab.
%
%  Although Matlab uses this representation, when we make graphs and
%  images we put the center of the image at the center -- of course -- and
%  we also put the DC value of the OTF in the middle of the image.  Hence,
%  when we return the frequency support or the spatial support we create
%  values for frequencies that run from negative to positive.  Similarly,
%  when we compute the spatial support we create spatial samples that run
%  below and above zero.
%
% Example:
%   oi = oiCreate; optics = oiGet(oi,'optics'); 
%   oi = oiSet(oi,'wave',400:10:700);
%
%   NA = opticsGet(optics,'na');   % Numerical aperture
%   rt = opticsGet(optics,'ray Trace');
%   psf = opticsGet(optics,'rtPSF',500);     % Shift-variant ray trace
%   psf = opticsGet(optics,'psf Data',600);  % Shift invariant data
%   vcNewGraphWin; mesh(sSupport(:,:,1),sSupport(:,:,2),psf);
%
%   otf = opticsGet(optics,'otf data',oi, 'mm',450); 
%   vcNewGraphWin; mesh(fftshift(abs(otf)));
%         
%   otfAll = opticsGet(optics,'otf data',oi);
%
% This OTF support does not work for ray trace optics (yet).
%   otfSupport = oiGet(oi,'fsupport','mm');  % Cycles/mm
%   vcNewGraphWin; mesh(otfSupport(:,:,1),otfSupport(:,:,2),fftshift(abs(otf)))
%
% Optics parameters
%
%   '*' means that you can use the syntax opticsGet(optics,'parm','unit'),
%   such as  opticsGet(optics,'focalLength','mm')
%
%      {'name'}    - name for these optics
%      {'type'}    - always 'optics'
%      {'model'}   - Type of optics computation,
%                    diffractionLimited, rayTrace, or shiftInvariant.
%      {'fnumber'}            - f#, ratio of focal length to aperture,
%                               a dimensionless quantity.
%      {'effective fnumber'}   - effective f-number
%      {'focal length'}        - focal length (M)
%      {'power'}               - optical power in diopters (1/f),units 1/M
%      {'image distance'}      - image distance from lensmaker's equation
%      {'image height'}*       - image height
%      {'imagew idth'}*        - image width
%          opticsGet(optics,'imagewidth',10,'mm')
%      {'image diagonal'}*     - image diagonal size
%      {'numerical aperture'}  - numerical aperture
%      {'aperture diameter'}*  - aperture diameter
%      {'aperture radius'}*    - aperture radius
%      {'aperture area'}*      - aperture area
%      {'magnification'}       - optical image magnification (<0 inverted)
%      {'pupil magnification'} -
%
% Off-axis methods and data
%      {'off axis method'}     - custom relative illumination method
%      {'cos4th method'}       - default cos4th method
%      {'cos4th data'}         - place to store cos4th data
%
% OTF information   - Used for shift-invariant calculations
%      {'otf data'}        - the optical transfer function data
%      {'otf size'}
%      {'otf fx'}          - column (fx) samples of OTF data
%      {'otf fy'}          - row (fy) samples of OTF data
%      {'otf support'}     - cell array, val{1:2}, of fy,fx samples
%      {'otf wave'}        - wavelength samples of the otf data
%      {'otf binwidth'}    - difference between wavelength samples
%      {'psf data'}        - psf data, calculated from the stored otfdata
%      {'psf spacing'}
%      {'psf support'}
%      {'incoherentcutoffspatialfrequency'}*    - Vector of incoherent cutoff freq
%                                                 for all wavelengths
%      {'maxincoherentcutoffspatialfrequency'}* - Largest incoherent cutoff
%
% Wavlength information
%      {'spectrum'}         - wavelength information
%        {'wavelength'}     - wavelength samples
%        {'nwave'}          - number of wavelength samples
%        {'binwidth'}       - spacing between the samples
%      {'transmittance'}    - wavelength transmission function
%
%  Ray Trace information - Used for non shift-invariant calculations
%      {'rtname'}        - name, may differ from file because of processing
%      {'raytrace'}      - structure of ray trace information
%      {'rtopticsprogram'}     - 'zemax' or 'code v'
%      {'rtlensfile'}          - Name of lens description file
%      {'rteffectivefnumber'}  - Effective fnumber
%      {'rtfnumber'}           - F-number
%      {'rtmagnification'}     - Magnification
%      {'rtreferencewavelength'}    - Design reference wavelength (nm)
%      {'rtobjectdistance'}*        - Design distance to object plane
%      {'rtfieldofview'}            - Diagonal field of view (deg)
%      {'rteffectivefocallength'}*  - Effective focal length
%      {'rtpsf'}               - structure containing psf information
%        {'rtpsfdata'}            - psf data
%                opticsGet(optics,'rtpsfdata')
%                opticsGet(optics,'rtpsfdata',fieldHeight,wavelength)
%        {'rtpsfsize'}            - (row,col) of the psf functions
%        {'rtpsfwavelength'}      - sample wavelengths of psf estimates
%        {'rtpsffieldheight'}*    - field heights for the psfs
%        {'rtpsfsamplespacing'}*  - sample spacing within the psfs
%        {'rtpsfsupport'}*        - spatial position (2D) of the psf functions
%        {'rtpsfsupportrow'}*     - spatial position of row samples
%        {'rtpsfsupportcol'}*     - spatial position of col samples
%        {'rtotfdata'}            - OTF derived from PSF ray trace data  *** (NYI)
%      {'rtrelillum'}       - structure of relative illumination information
%        {'rtrifunction'}       - Relative illumination function
%        {'rtriwavelength'}     - Wavelength samples (nm)
%        {'rtrifieldheight'}*   - Field heigh values
%      {'rtgeometry'}       - structure of geometric distortion information
%        {'rtgeomfunction'}         - Geometric distortion function
%               opticsGet(optics,'rtgeomfunction',[],'mm')
%               opticsGet(optics,'rtgeomfunction',500)
%        {'rtgeomwavelength'}       - Wavelength samples (nm)
%        {'rtgeomfieldheight'}*     - Field height samples
%        {'rtgeommaxfieldheight'}*  - Maximum field height sample
%
% Computational parameters
%       {'rtComputeSpacing'}*      - Sample spacing for PSF calculation
%
% Copyright ImagEval Consultants, LLC, 2005.


% Programming TODO:
%   Many of the rt spatial variables are stored in mm by the rtImportData
% function.  So, we are always dividing them by 1000 for return in meters.
% We should probably just store them in meters properly inside of
% rtImportData.
%
% The OTF support functions are in bad shape.  Fix.
%
val = [];

if ~exist('optics','var') || isempty(optics), 
    error('No optics specified.'); 
end
if ~exist('parm','var')   || isempty(parm), 
    error('No parameter specified.');
end

% We return different parameters depending on whether the user has a
% shift-invariant lens model (e.g., diffraction-limited) or a general ray
% trace model.
rt = 0;
if checkfields(optics,'rayTrace') && ~isempty(optics.rayTrace)
    % If there are ray trace data, and the current model is ray trace, 
    % set rt to 1.
    if strcmpi(optics.model,'raytrace'), rt = 1; end
end

parm = ieParamFormat(parm);
switch parm
    case 'name'
        val = optics.name;
    case 'type'
        val = optics.type;  % Should always be 'optics'
        
    case {'fnumber','f#'}
        % This is the f# assuming an object is infinitely far away.
        if rt
            if checkfields(optics,'rayTrace','fNumber'), val = optics.rayTrace.fNumber; end
        else val = optics.fNumber;
        end
    case {'model','opticsmodel'}
        if checkfields(optics,'model'), val = optics.model;
        else val = 'diffractionLimited';
        end
        
        % The user can set 'Diffraction limited' but have
        % 'diffractionlimited' returned.
        val = ieParamFormat(val);
        
    case {'effectivefnumber','efffnumber','efff#'}
        % The f# if the object is not at infinity.
        if rt
            if checkfields(optics,'rayTrace','effectiveFNumber'),
                val = optics.rayTrace.effectiveFNumber;
            end
        else
            val = opticsGet(optics,'fNumber')*(1 - opticsGet(optics,'mag'));
        end
    case {'focallength','flength'}
        % opticsGet(optics,'flength',units);
        if rt,  val = opticsGet(optics, 'RTeffectiveFocalLength');
        elseif strcmpi(opticsGet(optics,'model'),'skip')
 
            % If you choose 'skip' because you want to treat the
            % optics/lens as a pinhole, you must have a scene and in that
            % case we use the proper distance (half the scene distance). 
            % When you are just skipping to save time, you may not have a
            % scene. In that case, use the optics focal length.
            scene = vcGetObject('scene');
            if isempty(scene), val = optics.focalLength;
            else               val = sceneGet(scene,'distance')/2;
            end
            
        else   val = optics.focalLength;
        end
        if ~isempty(varargin)
            val = val*ieUnitScaleFactor(varargin{1});
        end

    case {'power','diopters'}
        % opticsGet(optics,'power','mm')
        % Diopters (1/m).
        % Sometimes we ask for a different unit when we don't want
        % diopters, but we want dist per deg in another unit.
        units = 'm';
        if ~isempty(varargin), units = varargin{1}; end
        val = 1/opticsGet(optics,'focallength',units);

    case {'imagedistance','focalplane','focalplanedistance'}
        % Lensmaker's equation calculation of focal plane distance from
        % the center of the lens.  If no source distance is provided, we
        % assume infinite distance
        %
        % opticsGet(optics,'focalplane',sDist); -- sDist is sourceDistance
        % opticsGet(optics,'focalplanedistance');  -- Infinite source dist

        % No need to check rt because focalLength parameter does
        % What about 'skip' case?  Should 'focalLength' call set this to
        % 1/2 the scene distance, which preserves the geometry?
        fL = opticsGet(optics,'focal Length');
        if isempty(varargin), sDist = Inf;
        else                  sDist = varargin{1};
        end
        val = 1 / (1/fL - 1/sDist);   % Lens maker equation

        % These all need the scene field of view.  They can be computed
        % from the geometry of the image distance and FOV.
    case {'imageheight'}
        % opticsGet(optics,'imageheight',fov) -- field of view in degrees
        if isempty(varargin), disp('fov required.'); return;
        else
            fov = varargin{1};
            imageDistance = opticsGet(optics,'focal plane distance');
            val = 2*imageDistance*tand(fov/2);
            if length(varargin) < 2, return;
            else val = ieUnitScaleFactor(varargin{2})*val;
            end
        end

    case {'imagewidth'}
        % opticsGet(optics,'imagewidth',fov) -- fov in degrees
        % opticsGet(optics,'imagewidth',fov,'mm') -- fov in degrees, output
        % units in mm
        if isempty(varargin), return;
        else
            fov = varargin{1};
            imageDistance = opticsGet(optics,'focalplanedistance');
            val = 2*imageDistance*tand(fov/2);
            if length(varargin) < 2, return;
            else val = ieUnitScaleFactor(varargin{2})*val;
            end
        end

    case {'imagediagonal','diagonal'}
        % opticsGet(optics,'imagediagonal',fov)  -- fov in degrees
        if isempty(varargin), return;
        else
            fov = varargin{1};
            h = opticsGet(optics,'imageheight',fov,'m');
            w = opticsGet(optics,'imagewidth',fov,'m');
            val = sqrt(h^2 + w^2);
        end
        if length(varargin) < 2, return;
        else val = ieUnitScaleFactor(varargin{2})*val;
        end

    case {'na','numericalaperture'}
        % Should this be a call to effective f#?
        val=1/(2*opticsGet(optics,'fnumber'));
    case {'aperturediameter','diameter','pupildiameter'}
        %These already check the rt condition, so no need to do it again
        val = opticsGet(optics,'focalLength')/opticsGet(optics,'fnumber');
        if ~isempty(varargin), val = ieUnitScaleFactor(varargin{1})*val; end
    case {'apertureradius','radius','pupilradius'}
        val = opticsGet(optics,'diameter')/2;
        if ~isempty(varargin), val = ieUnitScaleFactor(varargin{1})*val; end
    case {'aperturearea','pupilarea'}
        val = pi*((opticsGet(optics,'focalLength')/opticsGet(optics,'fnumber'))*(1/2))^2;
        if ~isempty(varargin), val = ieUnitScaleFactor(varargin{1})^2*val; end

    case {'magnification','mag'}
        % If the ray trace magnification is present, use that.  Otherwise
        % make a magnification estimate using the source distance and focal
        % length (via lensmaker equation).
        % opticsGet(optics,'mag',sDist) -- specify source distance
        % opticsGet(optics,'mag')       -- current source distance
        if rt
            val = opticsGet(optics,'rtmagnification');
            return;
        end
        
        % Skip model has unit magnification, but still negative
        if strcmpi(opticsGet(optics,'model'),'skip')
            val = -1;
            return;
        elseif length(varargin) < 1
            scene = vcGetObject('scene');
            if ~isempty(scene),  sDist = sceneGet(scene,'distance');
            else
                % What should we do here?  No scene and no distance
                % specified.  Do we make the magnification 0?  Seems OK.
                val = 0;  return;
            end
        else
            sDist = varargin{1};
        end
        % Distance to the image divided by distance to the object (source)
        % http://en.wikipedia.org/wiki/Magnification
        % Also f / (f - distObj)
        % If you want the mag to be -2, then 
        % distObj -2*f + 2*distObj = f, 2*distObj = 3*f, distObj = 3f/2
        % and in general distObj = (-M+1)f/(-M) * f
        val = -(opticsGet(optics,'focalPlaneDistance',sDist))/sDist;

    case {'pupilmagnification','pupmag'}
        % Pupil magnification is the ratio of the exit pupil to the input
        % pupil (diameters)
        if rt
            % Needs to be checked more thoroughly!  Based on PC.
            efn = opticsGet(optics,'rteffectiveFnumber');
            fn = opticsGet(optics,'fnumber');
            val = (1/(1 - (efn/fn)))*opticsGet(optics,'rtmagnification');
        else val = 1;
        end

    case {'wave','wavelength','wavelengthvalues'}
        % Wavelength related
        % nanometers is default
        % opticsGet(optics,'wavelength',unit)
        %
        if checkfields(optics,'spectrum','wave'), val = optics.spectrum.wave; end
        if isempty(val), scene = vcGetObject('scene'); val = sceneGet(scene,'wavelength'); end
        if isempty(val), val = 400:10:700; end
        
        if ~isempty(varargin)
            s = ieUnitScaleFactor(varargin{1})/ieUnitScaleFactor('nm');
            val = val*s;
        end
    case {'nwave','numberofwavelengthsamples'}
        if checkfields(optics,'spectrum','wave'), val = length(optics.spectrum.wave);
        end
    case {'binwidth','wavelengthbinwidth'}
        % Nanometers
        wave = opticsGet(optics,'wave');
        if length(wave) > 1, val = wave(2) - wave(1);
        else val = 1;
        end

    case {'transmittance','wavelengthtransmittance'}
        % [0,1]
        if checkfields(optics,'transmittance'), val = optics.transmittance;
        else val = ones(1,opticsGet(optics,'nwave')); end

        % ----- Diffraction limited parameters
    case {'dlfsupport','dlfsupportmatrix'}
        % Two different return formats.  Either
        %  val{1} and val{2} as vectors, or
        %  val  = fSupport(:,:,:);
        % opticsGet(optics,'dl fsupport',wave,unit,nSamp)
        % opticsGet(optics,'dl fsupport matrix',wave,unit,nSamp)
        % 
        % Diffraction limited frequency support at a wavelength (i.e.,
        % support out to the incoherent cutoff frequency).  This can be
        % used for plotting, for example.

        if length(varargin) < 1, error('Must specify wavelength'); else thisWave = varargin{1}; end
        if length(varargin) < 2, units = 'mm'; else units = varargin{2}; end
        if length(varargin) < 3, nSamp = 30; else nSamp = varargin{3}; end

        % Sometimes the optics wavelength hasn't been defined because, say,
        % we haven't run through a scene.  So we trap that case here.
        waveList = opticsGet(optics,'wavelength');
        idx  = ieFindWaveIndex(waveList,thisWave);
        inCutoff = opticsGet(optics,'inCutoff',units); 
        inCutoff = inCutoff(idx);
        
        fSamp = (-nSamp:(nSamp-1))/nSamp;
        val{1} = fSamp*inCutoff;
        val{2} = fSamp*inCutoff;

        % Alternative return format
        if strfind(parm,'matrix')
            [valMatrix(:,:,1),valMatrix(:,:,2)] = meshgrid(val{1},val{2});
            val = valMatrix;
        end
        
    case {'incoherentcutoffspatialfrequency','incutfreq','incutoff'}
        % cycles/distance
        % Cutoff spatial frequency for a diffraction limited lens.  See
        % formulae in dlCore.m
        apertureDiameter = opticsGet(optics,'aperturediameter');
        imageDistance    = opticsGet(optics,'focalplanedistance');
        wavelength       = opticsGet(optics,'wavelength','meters');

        % Sometimes the optics wavelength have not been assigned because
        % there is no scene and no oiCompute has been run.  So, we can just
        % choose a sample set.
        if isempty(wavelength), wavelength = (400:10:700)*10^-9; end

        % See dlCore.m for a description of the formula.  We divide by the
        % scale factor, instead of multiplying, because these are
        % frequencies (1/m), not distances.
        val = (apertureDiameter / imageDistance) ./ wavelength;
        if ~isempty(varargin), val = val/ieUnitScaleFactor(varargin{1}); end

    case {'maxincoherentcutoffspatialfrequency','maxincutfreq','maxincutoff'}
        % opticsGet(optics,'maxincutoff','m')
        % opticsGet(optics,'maxincutoff')
        if isempty(varargin), val = max(opticsGet(optics,'incutoff'));
        else val = max(opticsGet(optics,'incutoff',varargin{1}));
        end

        % -------   OTF information and specifications.  This case is used for
        % shift-invariant calculations.  The ray trace structures below are
        % used for non shift-invariant cases derived from Zemax or Code V data
        % sets.

    case {'otf','otfdata','opticaltransferfunction'}
        % You can ask for a particular wavelength with the syntax
        %    opticsGet(optics,'otfData',oi, spatialUnits, wave)
        %
        % OTF values can be complex. They are related to the PSF data by
        %    OTF(:,:,wave) = fft2(psf(:,:,wave))
        % For example, the PSF at 450 nm can be obtained via
        %    mesh(abs(fft2(opticsGet(optics,'otfdata',450))))
        %
        % We are having some issues on this point for shift-invariant and
        % diffraction limited models.  Apparently there is a problem with
        % fftshift???
       
        opticsModel = opticsGet(optics,'model');
        thisWave = [];
        switch lower(opticsModel)
            case 'diffractionlimited'
                % For diffraction limited case, the call must be 
                % otf = opticsGet(optics,'otf data',oi, fSupport, [wave]);
                units = 'mm';
                if isempty(varargin)
                    error('format is ... oi, spatialUnits,wave');
                else oi = varargin{1};
                end
                
                if length(varargin) > 1, units = varargin{2};
                end
                if length(varargin) > 2, thisWave = varargin{3}; end
                
                fSupport = oiGet(oi,'fSupport',units);   % 'cycles/mm'
                wavelength = oiGet(oi,'wave');
                OTF = dlMTF(oi,fSupport,wavelength,'millimeters');

            case 'shiftinvariant'
                % opticsGet(optics,'otf data',[wave]);
                if checkfields(optics,'OTF','OTF')
                    OTF = optics.OTF.OTF;
                    if ~isempty(varargin), thisWave = varargin{1}; end
                else OTF = []; 
                end
                
            case 'raytrace'
                error('opticsGet(optics,''OTF'') not supported for ray trace');
                
            otherwise
                error('OTFData not implemented for %s model',opticsModel);
        end

        % Wavelength is asked for
        if ~isempty(thisWave)
            [idx1,idx2] = ieWave2Index(opticsGet(optics,'otfWave'),thisWave);
            if idx1 == idx2
                val = OTF(:,:,idx1);
            else
                wave = opticsGet(optics,'otfwave');
                w   = 1 - ((varargin{1} - wave(idx1))/(wave(idx2) - wave(idx1)));
                val = (w*OTF(:,:,idx1) + (1-w)*OTF(:,:,idx2));
                % wave(idx1),  varargin{1}, wave(idx2), w
            end
        else
            % Returns the entire data set
            val = OTF;
        end

    case {'psf'}
        % psf = opticsGet(optics,'psf',thisWave,units,nSamp,freqOverSample);
        % Pointspread at a particular wavelength, spatial sampling in some
        % units ('um','mm', etc.) some number of spatial samples, and some
        % amount of oversampling on the frequency calculation to make the
        % curve smooth.
        %
        if isempty(varargin), error('You must specify wavelength'); 
        else   thisWave = varargin{1};
        end
        if length(varargin) < 2, units = 'um'; else units = varargin{2}; end
        if length(varargin) < 3, nSamp = 100; else nSamp = varargin{3}; end
        if length(varargin) < 4, oSample = 1; else oSample = varargin{4}; end
        
        fSupport = opticsGet(optics,'dl fsupport matrix',thisWave,units,nSamp);
        fSupport = fSupport*oSample;        

        %  Oversample the frequency to get a smoother PSF image.
        %  You can specify the factor for oversampling in the
        %  calling arguments.
        otf = dlMTF(optics,fSupport,thisWave,units);
        val = fftshift(ifft2(otf));
        
    case {'degreesperdistance','degperdist'}
        % opticsGet(optics,'deg per dist','mm')
        % We use this constant to convert from the input spatial frequency units
        % (cycles/deg) to cycles/meter needed for the Hopkins eye. We need to
        % calculate this value from the optics data passed in.
        %
        % Given D0, the focal plane is 1/D0 meters from the lens.  Call this the
        % adjacent edge of the right triangle from the image plane to the lens.
        %
        % 1 deg of visual angle is
        %   tan(opp/(1/D0)) = 1             (deg)
        %   opp/(1/D0) = atand(1)           (1/rad)
        %   opp = atand(1)*(1/D0)           (1/rad * meter)
        %   1/opp = 1/ (atand(1)*(1/D0))    (rad/meter)
        %
        % The conversion is: (cycles/rad) * (rad/meter) = cycles/meter
        units = 'm'; 
        if ~isempty(varargin), units = varargin{1}; end
        D0 = opticsGet(optics,'power',units); 
        val = 1/(1/D0 * atand(1));  
        
    case {'distperdeg','distanceperdegree'}
        units = 'm';
        if ~isempty(varargin), units = varargin{1}; end
        val = 1/opticsGet(optics,'deg per dist',units);
        
        
        %------------------------
        % I think we should delete the stuff below and the OTF.fx entry,
        % using only the oiFrequencySupport call.  This will require
        % debugging, particularly in the script s_TestSI and perhaps other
        % places.
    case {'otfsupport'}
        % val = opticsGet(optics,'otf support','mm');
        %
        % Row and col (Y,X) spatial frequency range for the OTF data
        % [Y,X] meshgrid(val{1},val{2}) produces the matrices for surface
        % plotting.  
        % Frequency is stored in non-standard units of cycles/mm. This will
        % be annoying to fix some day, sigh.
        units = 'mm';
        if ~isempty(varargin), units = varargin{1}; end
        val{1} = opticsGet(optics,'otf fy',units);
        val{2} = opticsGet(optics,'otf fx',units);

    case {'otffx'}
        % cycles/mm!!! Non-standard unit. Must fix up some day.
        if checkfields(optics,'OTF','fx'), val = optics.OTF.fx; end
        % Transform into other units if required
        if ~isempty(varargin)
            unit = ieParamFormat(varargin{1});
            if strcmp(unit, 'cycles/deg') || strcmp(unit, 'cyclesperdeg')
                val = val*tand(1)*opticsGet(optics, 'focal length', 'mm');
            else
                val = (val*10^3)/ieUnitScaleFactor(unit);
            end
        end
    case {'otffy'}
        % cycles/mm!!! Non-standard unit. Must fix up some day.
        if checkfields(optics,'OTF','fy'), val= optics.OTF.fy; end
        % Put into meters and then apply scale factor
        if ~isempty(varargin), 
            unit = ieParamFormat(varargin{1});
            if strcmp(unit, 'cycles/deg') || strcmp(unit, 'cyclesperdeg')
                val = val*tand(1)*opticsGet(optics, 'focal length', 'mm');
            else
                val = (val*10^3)/ieUnitScaleFactor(unit);
            end
        end
    case {'otfsize'}
        % Row and col samples
        if checkfields(optics,'OTF','OTF'),
            tmp = size(optics.OTF.OTF); val = tmp(1:2);
        end

    case {'otfwave'}
        % opticsGet(optics,'otf wave','nm');
        % nm is the default.
        % This should probably go away and we should only wave 'wave'.
        if checkfields(optics,'OTF','wave'), val = optics.OTF.wave; 
        else val = opticsGet(optics,'wave');
        end
        if ~isempty(varargin)
            units = varargin{1};
            val = val*1e-9*ieUnitScaleFactor(units);
        end
        
    case {'otfbinwidth'}
        otfWave = opticsGet(optics,'otfWave');
        if length(otfWave)>1, val = otfWave(2) - otfWave(1);
        else val = 1;
        end
    case {'psfdata'}
        % To return the psf at 500 nm use
        %    psf = opticsGet(optics,'psfData',500);
        %    mesh(psf);
        % The 
        oModel = opticsGet(optics,'model');
        switch lower(oModel)
            case 'diffractionlimited'
                % opticsGet(optics,'psf Data',500,'um');
                if length(varargin) < 1
                    % Person wanted all wavelengths
                    thisWave = opticsGet(optics,'wave');
                else thisWave = varargin{1};
                end
                if length(varargin) < 2, units = 'um';
                else units = varargin{2};
                end
                
                nSamp = 100;   % Number of frequency steps from 0 to incoherent cutoff
                dlF = opticsGet(optics,'dlFSupport',thisWave(1),units,nSamp);
                [fSupport(:,:,1),fSupport(:,:,2)] = meshgrid(dlF{1},dlF{2});
                % This is a way to calculate the spatial support for the
                % psf
                % See 'psfsupport' below, also.  These should be integrated
                % and coordinated with fSupport.
                % If we would like the spatial support to be smaller
                % and finer, we should scale fSupport*4, above.
                %   samp = (-nSamp:(nSamp-1));
                %   [X,Y] = meshgrid(samp,samp);
                %   deltaSpace = 1/(2*max(fSupport(:)));
                %   sSupport(:,:,1) = X*deltaSpace;
                %   sSupport(:,:,2) = Y*deltaSpace;

                % Diffraction limited OTF
                otf = dlMTF(optics,fSupport,thisWave,units);
                if length(thisWave) == 1
                    val = fftshift(ifft2(otf));
                    val = abs(val);
                else
                    val = zeros(size(otf));
                    for ii=1:length(thisWave)
                        val(:,:,ii) = fftshift(ifft2(otf(:,:,ii)));
                    end
                end

            case 'shiftinvariant'
                if checkfields(optics,'OTF','OTF')
                    otfWave = opticsGet(optics,'otfWave');
                    if ~isempty(varargin)
                        % Just at the interpolated wavelength
                        val = opticsGet(optics,'otfData',varargin{1});
                        val = fftshift(ifft2(val));
                        % mesh(val)
                    else
                        % All of them
                        val = zeros(size(optics.OTF.OTF));
                        for ii=1:length(otfWave)
                            val(:,:,ii) = fftshift(ifft2(optics.OTF.OTF(:,:,ii)));
                        end
                    end
                    if ~isreal(val)
                        warning('ISET:complexpsf','complex psf');
                        val = abs(val);
                    end
                else
                    warning('ISET:otfdata','No OTF data stored in optics.')
                end
        end

    case {'psfspacing'}
        % opticsGet(optics,'psf spacing',unit)
        % Sample spacing of the psf points
        %
        % Warning:  We are assuming that fx and fy have the same peak
        % spatial frequency and spatial sampling.
        if length(varargin) >= 1, units = varargin{1}; 
        else units = 'mm'; end
        fx = opticsGet(optics,'otf fx',units); 
        peakF = max(fx(:));
        if isempty(fx), error('No otffx calculated yet. Fix me.'); end

        % Peak frequency in cycles/meter.  1/peakF is meters.  We have two
        % samples in that distance, so the sample spacing is half that
        % distance.
        val = 1/(2*peakF);
    case {'psfsupport'}
        % opticsGet(optics,'psf support',unit)
        % Returns mesh grid of X and Y values.  Used for mesh plotting
        % often.
        % X/Y could be mixed up in 1 and 2.
        % This should be replaced by the code in plotOI/OTF psf550 case.
        
        if length(varargin) >= 1, units = varargin{1}; 
        else units = 'mm'; end
        
        sz = opticsGet(optics,'otf size');
        if isempty(sz), error('No optical image data'); end
        
        x = (0:(sz(1)-1))*opticsGet(optics,'psfspacing',units);
        x = x - mean(x);
        [X,Y] = meshgrid(x,x);
        val{1} = X; val{2} = Y;

        %----------- Relative illumination (off-axis) specifications
    case {'offaxis','offaxismethod','relativeilluminationtype'}
        % This is the method used to compute relative illumination. It can
        % be 'Skip','cos4th','codeV', or 'Zemax'.
        val = optics.offaxis;
    case {'cos4thmethod','cos4thfunction'}
        % Most people use cos4th as an offaxis method. In that case, the
        % function that implements cos4th can be stored here.  I suspect
        % this extra step is not needed. We have a cos4th function and we use
        % that without allowing some other implementation.  It is here only
        % for some old backwards compatibility.
        %
        % Do not run cos4th when you are using the Code V (or probably
        % Zemax methods.  These calculations include the cos4th mechanisms
        % in the lens calculations.
        if checkfields(optics,'cos4th','function'), val = optics.cos4th.function; end
    case {'cos4th','cos4thdata','cos4thvalue'}
        % Numerical values.  Should change field to data from value.  I
        % don't think this is ever used, is it?
        if checkfields(optics,'cos4th','value'), val = optics.cos4th.value; end

        % ---------------  Ray Trace information.
        % The ray trace computations differ from those above because they
        % are not shift-invariant.  When we use a custom PSF/OTF that is
        % shift invariant, we still store the information in the main
        % optics code region in the OTF structure.
    case {'raytrace','rt',}
        if checkfields(optics,'rayTrace'), val = optics.rayTrace; end
    case {'rtname'}
        if checkfields(optics,'rayTrace','name'), val = optics.rayTrace.name; end
    case {'opticsprogram','rtopticsprogram'}
        if checkfields(optics,'rayTrace','program'), val = optics.rayTrace.program; end
    case {'lensfile','rtlensfile'}
        if checkfields(optics,'rayTrace','lensFile'), val = optics.rayTrace.lensFile;end

    case {'rteffectivefnumber','rtefff#'}
        if checkfields(optics,'rayTrace','effectiveFNumber'), val = optics.rayTrace.effectiveFNumber;end
    case {'rtfnumber'}
        if checkfields(optics,'rayTrace','fNumber'), val = optics.rayTrace.fNumber;end
    case {'rtmagnification','rtmag'}
        if checkfields(optics,'rayTrace','mag'), val = optics.rayTrace.mag; end
    case {'rtreferencewavelength','rtrefwave'}
        if checkfields(optics,'rayTrace','referenceWavelength'), val = optics.rayTrace.referenceWavelength;end
    case {'rtobjectdistance','rtobjdist','rtrefobjdist','rtreferenceobjectdistance'}
        % TODO:  These are stored in mm, I believe.  Could change to m
        if checkfields(optics,'rayTrace','objectDistance'),
            val = optics.rayTrace.objectDistance/1000;
        end
        if ~isempty(varargin), val = val*ieUnitScaleFactor(varargin{1}); end

    case {'rtfieldofview','rtfov','rtmaximumfieldofview','rtdiagonalfov', 'rthorizontalfov'}
        % Maximum diagonal field of view for the ray trace calculation (not
        % the computed image).
        %
        % The stored value is half the maximum diagonal.  The max
        % horizontal is the same.  The FOV (whatever direction) is as far
        % as the PSF is computed by Zemax. It doesn't matter whether it is
        % measured along the diagonal or the horizontal.
        %
        % Until Aug. 13 2011, however, we stored sqrt(2)*fov, rather than
        % the fov.  This was an error at ImageVal.
        if checkfields(optics,'rayTrace','maxfov'), val = optics.rayTrace.maxfov; end
    case {'rteffectivefocallength','rtefl','rteffectivefl'}
        % TODO:  These are stored in mm, I believe.  Should change to m
        if checkfields(optics,'rayTrace','effectiveFocalLength'),
            val = optics.rayTrace.effectiveFocalLength/1000;
        end
        if ~isempty(varargin), val = val*ieUnitScaleFactor(varargin{1}); end

    case {'rtpsf'}
        if checkfields(optics,'rayTrace','psf'), val = optics.rayTrace.psf; end
    case {'rtpsffunction','rtpsfdata'}
        % Return the psf - either the whole thing or a selected psf
        % Data are stored as 4D (row,col,fieldHeight,wavelength) images
        if checkfields(optics,'rayTrace','psf','function'),
            if ~isempty(varargin)
                % Return the psf at a particular field height and wavelength
                % The units of the field height are meters and nanometers
                % for wavelength
                % psf = opticsGet(optics,'rtpsfdata',fieldHeight,wavelength);
                % Delete this warning after January 2009 (warning commented
                % out April 2011)
                % if varargin{1} > .1, warndlg('Suspiciously large field height (> 0.1m)'); end
                fhIdx   = ieFieldHeight2Index(opticsGet(optics,'rtPSFfieldHeight'),varargin{1});
                waveIdx = ieWave2Index(opticsGet(optics,'rtpsfwavelength'),varargin{2});
                val = optics.rayTrace.psf.function(:,:,fhIdx,waveIdx);
             else
                % Return the entire psf data
                % psfFunction = opticsGet(optics,'rtpsfdata');
                val = optics.rayTrace.psf.function;
            end
        end
    case {'rtpsfsize','rtpsfdimensions'}
        % psfSize = opticsGet(optics,'rtPsfSize')
        if checkfields(optics,'rayTrace','psf','function')
            % All 4 dimensions, in the order row,col,fieldnum,wavelength
            val = size(optics.rayTrace.psf.function);
        end
    case {'rtpsfwavelength'}
        if checkfields(optics,'rayTrace','psf','wavelength'), val = optics.rayTrace.psf.wavelength; end
    case {'rtpsffieldheight'}
        % opticsGet(optics,'rt psf field height','um')
        % Stored in mm. Returned in the requested units.
        if checkfields(optics,'rayTrace','psf','fieldHeight'),
            val = optics.rayTrace.psf.fieldHeight/1000;
        end
        if ~isempty(varargin), val = val*ieUnitScaleFactor(varargin{1}); end

    case {'rtpsfsamplespacing','rtpsfspacing'}
        % opticsGet(optics,'rtPsfSpacing','um')
        if checkfields(optics,'rayTrace','psf','sampleSpacing')
            % The 1000 is necessary because it is stored in mm
            val = optics.rayTrace.psf.sampleSpacing/1000;
        end 
        if ~isempty(varargin), val = val*ieUnitScaleFactor(varargin{1}); end
    case {'rtsupport','rtpsfsupport'}
        % Return the (x,y) positions of the PSF samples.
        % This list always contains a sample at 0
        % For a 2^n (e.g., 64,64) grid, Code V puts the zero point at 2^n - 1 (e.g., 33,33)
        % s = opticsGet(optics,'rtPsfSupport','um');
        psfSize = opticsGet(optics,'rtPSFSize');
        if isempty(varargin), units = 'm'; else units = varargin{1}; end
        psfSpacing = opticsGet(optics,'rtPsfSpacing',units);
        colPsfPos = (((-psfSize(2)/2)+1):(psfSize(2)/2))*psfSpacing(2);
        rowPsfPos = (((-psfSize(1)/2)+1):(psfSize(1)/2))*psfSpacing(1);
        [xPSF,yPSF] = meshgrid(colPsfPos,rowPsfPos);
        val(:,:,1) = xPSF;
        val(:,:,2) = yPSF;
    case {'rtpsfsupportrow','rtpsfsupporty'}
        psfSize = opticsGet(optics,'rtPSFSize');
        if isempty(varargin), units = 'm'; else units = varargin{1}; end
        psfSpacing = opticsGet(optics,'rtPsfSpacing',units);
        val = (((-psfSize(1)/2)+1):(psfSize(1)/2))*psfSpacing(1);
        % Useful to be a column for interpolation.  Maybe they should
        % always be columns?
        val = val(:);

    case {'rtpsfsupportcol','rtpsfsupportx'}
        psfSize = opticsGet(optics,'rtPSFSize');
        if isempty(varargin), units = 'm'; else units = varargin{1}; end
        psfSpacing = opticsGet(optics,'rtPsfSpacing',units);
        val = (((-psfSize(2)/2)+1):(psfSize(2)/2))*psfSpacing(2);
        % Useful to be a row for interpolation.  But maybe it should always
        % be a column?
        val = val(:)';
    case {'rtfreqsupportcol','rtfreqsupportx'}
        % Calculate the frequency support across the column dimension
        % opticsGet(optics,'rtFreqSupportX','mm')
        if isempty(varargin), units = 'm'; else units = varargin{1}; end
        psfSpacing = opticsGet(optics,'rtPsfSpacing',units);
        sz = opticsGet(optics,'rtPsfSize',units);
        val = (((-sz(2)/2)+1):(sz(2)/2))*(1/(sz(2)*psfSpacing(2)));
        val = val(:)';
    case {'rtfreqsupportrow','rtfreqsupporty'}
        % Calculate the frequency support across the column dimension
        % opticsGet(optics,'rtFreqSupportX','mm')
        if isempty(varargin), units = 'm'; else units = varargin{1}; end
        psfSpacing = opticsGet(optics,'rtPsfSpacing',units);
        sz = opticsGet(optics,'rtPsfSize',units);
        val = (((-sz(1)/2)+1):(sz(1)/2))*(1/(sz(1)*psfSpacing(1)));
        val = val(:);
    case {'rtfreqsupport'}
        % val = opticsGet(optics,'rtFreqSupport','mm');
        if isempty(varargin), units = 'm'; else units = varargin{1}; end
        val{1} = opticsGet(optics,'rtFreqSupportX','mm');
        val{2} = opticsGet(optics,'rtFreqSupportY','mm');
    case {'rtrelillum'}
        % The sample spacing on this is given below in rtrifieldheight
        if checkfields(optics,'rayTrace','relIllum'), val = optics.rayTrace.relIllum; end
    case {'rtrifunction','rtrelativeilluminationfunction','rtrelillumfunction'}
        if checkfields(optics,'rayTrace','relIllum','function'),
            val = optics.rayTrace.relIllum.function;
        end
    case {'rtriwavelength','rtrelativeilluminationwavelength'}
        if checkfields(optics,'rayTrace','relIllum','wavelength'),
            val = optics.rayTrace.relIllum.wavelength;
        end
    case {'rtrifieldheight','rtrelativeilluminationfieldheight'}
        % TODO:  These are stored in mm, I believe.  Could change to m
        if checkfields(optics,'rayTrace','relIllum','fieldHeight'),
            val = optics.rayTrace.relIllum.fieldHeight/1000;
        end
        if ~isempty(varargin), val = val*ieUnitScaleFactor(varargin{1}); end

    case {'rtgeometry'}
        if checkfields(optics,'rayTrace','geometry'), val = optics.rayTrace.geometry; end
    case {'rtgeomfunction','rtgeometryfunction','rtdistortionfunction','rtgeomdistortion'}
        % opticsGet(optics,'rtDistortionFunction',wavelength,'units');
        if checkfields(optics,'rayTrace','geometry','function'),
            % opticsGet(optics,'rtGeomFunction',[],'mm')
            % opticsGet(.,.,500,'um')  % 500 nm returned in microns
            % opticsGet(.,.)           % Whole function in meters
            % opticsGet(optics,'rtgeomfunction',500) % Meters
            if isempty(varargin)
                % Return the whole function units are millimeters
                val = optics.rayTrace.geometry.function;
                return;
            else
                % Return values at a specific wavelength
                if ~isempty(varargin{1})
                    idx = ieWave2Index(opticsGet(optics,'rtgeomwavelength'),varargin{1});
                    val = optics.rayTrace.geometry.function(:,idx);
                else
                    val = optics.rayTrace.geometry.function;
                end
            end
   
            % Stored in millimeters. Convert to meters.
            val = val/1000;
            % If there is a second varargin, it specifieds the units.
            if length(varargin) == 2,
                val = val*ieUnitScaleFactor(varargin{2});
            end
        end

    case {'rtgeomwavelength','rtgeometrywavelength'}
        % The wavelength used for ray trace geometry distortions. 
        % The units is nanometers
        if checkfields(optics,'rayTrace','geometry','wavelength')
            val = optics.rayTrace.geometry.wavelength; 
        end
    case {'rtgeomfieldheight','rtgeometryfieldheight'}
        % val = opticsGet(optics,'rtGeomFieldHeight','mm');
        % These are stored in mm because of Zemax.  So we divide by 1000 to
        % put the value into meters and then convert to the user's
        % requested units.
        if checkfields(optics,'rayTrace','geometry','fieldHeight')
            % Convert from mm to meters
            val = optics.rayTrace.geometry.fieldHeight/1000;
        end
        if ~isempty(varargin), val = val*ieUnitScaleFactor(varargin{1}); end

    case {'rtgeommaxfieldheight','rtmaximumfieldheight','rtmaxfieldheight'}
        % val = opticsGet(optics,'rtGeomMaxFieldHeight','mm');
        % The maximum field height.
        fh = opticsGet(optics,'rtgeometryfieldheight');  % Returned in meters
        val = max(fh(:));
        if ~isempty(varargin), val = val*ieUnitScaleFactor(varargin{1}); end

    case {'rtcomputespacing'}
        % Spacing of the point spread function samples.
        % TODO: Is this really stored in meters, not millimeters like other
        % stuff?
        if checkfields(optics,'rayTrace','computation','psfSpacing'),
            val = optics.rayTrace.computation.psfSpacing;
            if ~isempty(varargin), val = val*ieUnitScaleFactor(varargin{1}); end
        end

    otherwise
        error('Unknown optics parameter.');

end

end