function oi = oiSet(oi,parm,val,varargin)
% Set ISET optical image parameter values
%
%    oi = oiSet(oi,parm,val,varargin)
%
% The parameters of an optical image are set through the calls to this
% routine.  Parameters of the optics, attached to this structure, can also
% be set via this call or by retrieving the optics structure and setting it
% directly (see below).
%
% The oi is the object; parm is the name of the parameter; val is the
% value of the parameter; varargin allows some additional parameters in
% certain cases.
%
% There is a corresponding oiGet routine.  Many fewer parameters are
% available for 'oiSet' than 'oiGet'. This is because many of the
% parameters returned by oiGet are derived from the few parameters that can
% be set.
%
% It is  possible to set to the values of the optics structure attached to
% the oi structure.  To do this, use the syntax
%
%   oiSet(oi,'optics param',val) where param is the optics parameter.
%
% This synatx replaces the older and more tedious style  
%    optics = oiGet(oi,'optics'); 
%    optics = opticsSet(optics,'param',val);
%    oi = oiSet(oi,'optics',optics);
%
% Examples:
%    oi = oiSet(oi,'name','myName')
%    oi = oiSet(oi,'filename','test')
%    oi = oiSet(oi,'optics',optics);
%    oi = oiSet(oi,'optics fnumber',2.8); oiGet(oi,'optics fnumber')
%
% N.B.: Because of the large size of the photon data (row,col,wavelength)
% and the high dynamic range, they are stored in a special compressed
% format.  They also permit the user to read and write individual
% wavelength planes. Data sent in and returned are always in double()
% format.
%
% When you write to 'photons', the compression fields used by cphotons are
% cleared. When reading and writing a waveband in compressed mode, it is
% assumed that the compression fields already exist.  We do not compress
% each individual waveband, though this would be possible (i.e., to have an
% array of min/max values for each waveband).
%
% After writing to the photons field, the illuminance and mean illuminance
% fields are set to empty.
%
% User-settable oi parameters
%
%      {'name'}
%      {'type'}
%      {'distance' }
%      {'horizontal field of view'}
%      {'magnification'}
%
%      {'data'}  - Irradiance information
%        {'cphotons'}   - Compressed photons; can be set one waveband at a
%                         time: oi = oiSet(oi,'cphotons',data,wavelength);
%
% Wavelength information
%      {'spectrum'}            - Spectrum structure
%        {'wavelength'}        - Wavelength samples
%
% Optics
%      {'optics'}  - Main optics structure
%      {'optics model'} - Optics computation
%         One of raytrace, diffractionlimited, or shiftinvariant 
%         Spaces and case variation is allowed, i.e.
%         oiSet(oi,'optics model','diffraction limited');
%         The proper data must be loaded to run oiCompute.
%
%      {'diffuser method'} - 'blur', 'birefringent' or 'skip'
%      {'diffuser blur'}   - FWHM blur amount (meters)
%
%      {'psfstruct'}        - Entire PSF structure (shift-variant)
%       {'sampled RT psf'}     - Precomputed shift-variant psfs
%       {'psf sample angles'}  - Vector of sample angle
%       {'psf image heights'}  - Vector of sampled image heights
%       {'rayTrace optics name'}  - Optics used to derive shift-variant psf
%
%      {'depth map'}         - Distance of original scene pixel (meters)
%
% Auxiliary
%      {'consistency'}       - Is the display consistent with data
%      {'gamma'}             - Display gamma in oiWindow
%
% Private variables used by ISET but not normally set by the user
%
%   Used for management of compressed photons
%      {'bit depth'}
%
%   Used to cache optical image illuminance
%      {'illuminance'}
%      {'mean illuminance'}
%
% Copyright ImagEval Consultants, LLC, 2003.

if ~exist('parm','var') || isempty(parm), error('Param must be defined.'); end
if ~exist('val','var'), error('Value field required.'); end

[oType,parm] = ieParameterOtype(parm);

% Handling optics setting via oiSet call.
if isequal(oType,'optics')
    if isempty(parm)
        % oi = oiSet(oi,'optics',optics);
        oi.optics = val;
        return;
    else
        % Allows multiple additional arguments
        if isempty(varargin), oi.optics = opticsSet(oi.optics,parm,val);
        elseif length(varargin) == 1
            oi.optics = opticsSet(oi.optics,parm,val,varargin{1});
        elseif length(varargin) == 2
            oi.optics = opticsSet(oi.optics,parm,val,varargin{1},varargin{2});
        elseif length(varargin) == 3
            oi.optics = opticsSet(oi.optics,parm,val,varargin{1},varargin{2},varargin{3});
        end
        return;
    end
elseif isempty(parm)
    error('oType %s. Empty param.\n',oType);
end


parm = ieParamFormat(parm);

switch parm

    case {'name','oiname'}
        oi.name = val;
    case 'type'
        oi.type = val;
    case {'filename'}
        % When the data are ready from a file, we save the file name.
        % Happens, perhaps, when reading multispectral image data.
        oi.filename = val;
    case {'consistency','computationalconsistency'}
        % When parameters are changed, the consistency flag on the optical
        % image changes.  This is irrelevant for the scene case.
        oi.consistency = val;
    case {'gamma'}
        % oiSet([],'gamma',0.6);
        % Should this be ieSessionSet('oi gamma',val)
        hObj = ieSessionGet('oi window ');
        hdl = ieSessionGet('oi window handle');
        eventdata = [];
        set(hdl.editGamma,'string',num2str(val));
        oiWindow('oiRefresh',hObj,eventdata,hdl);
        
    case {'distance' }
        % Positive for scenes, negative for optical images
        oi.distance = val;

    case {'wangular','widthangular','hfov','horizontalfieldofview','fov'}
        oi.wAngular = val;

    case 'magnification'
        % Optical images have other mags calculated from the optics.
        evalin('caller','mfilename')
        warndlg('Setting oi magnification.  Bad idea.')
        oi.magnification = val;

    case {'optics','opticsstructure'}
        oi.optics = val;

    case {'data','datastructure'}
        oi.data = val;

    case {'photons','cphotons', 'compressedphotons'}
        % oiSet(oi,'photons',val)
        % cphotons is obsolete and should not be used.  It will go away
        % before long.
        if ~(isa(val, 'double') || isa(val, 'single') || isa(val, 'gpuArray')),
            error('Photons must be type double / single / gpuArray');
        end
        
        % Perhaps this should go away, too.  We should just stick with
        % single and gpuArray throughout.
        bitDepth = oiGet(oi, 'bitDepth');
        if isempty(bitDepth), error('Compression parameters not set'); end
        if ~isempty(varargin)
            % varargin{1} - selected waveband
            idx = ieFindWaveIndex(oiGet(oi, 'wave'),varargin{1});
            idx = logical(idx);
        end
        
        switch bitDepth
            case 64 % Double
                if isempty(varargin)
                    oi.data.photons = val;
                else
                    oi.data.photons(:,:,idx) = val;
                end
            case 32 % Single
                if isempty(varargin)
                    oi.data.photons = single(val);
                else
                    oi.data.photons(:,:,idx) = single(val);
                end
            otherwise
                error('Unsupported bit depth %f',bitDepth);
        end

        % Clear out derivative luminance/illuminance computations
        oi = oiSet(oi,'illuminance', []);
        oi = oiSet(oi,'meanilluminance', []);

    case {'datamin','dmin'}
        % Only used by compressed photons. Not by user
        error('datamin and datamax are not used anymore');
        % oi.data.dmin = val;
    case {'datamax','dmax'}
        % Only used by compressed photons. Not by user
        error('datamin and datamax are not used anymore');
        % oi.data.dmax = val;
    case 'bitdepth'
        % Only used by compressed photons
        oi.data.bitDepth = val;

    case {'illuminance', 'illum'}
        % The value is stored for efficiency.
        oi.data.illuminance = val;

    case {'meanillum', 'meanilluminance'}
        oi.data.meanIll = val;

    case {'spectrum','wavespectrum','wavelengthspectrumstructure'}
        oi.spectrum  = val;
        %     case {'binwidth','wavelengthspacing'}
        %         oi.spectrum.binwidth = val;
    case {'wave','wavelength','wavelengthnanometers'}
        % We should probably check that val is a proper set of wavelength
        % values that make sense ... unique, evenly spaced, stuff like
        % that.
        oi.spectrum.wave = val;

        % Optical methods
    case {'opticsmodel'}
        % oi = oiSet(oi,'optics model', 'ray trace');
        % The optics model should be one of
        % raytrace, diffractionlimited, or shiftinvariant
        % Spacing and case variation is allowed.
        val = ieParamFormat(val);
        oi.optics.model = val;

        % Glass diffuser properties
    case {'diffusermethod'}
        % This determines calculation 
        % 0 - skip, 1 - gaussian blur, 2 - birefringent
        % We haven't set up the interface yet (12/2009)
        oi.diffuser.method = val;
    case {'diffuserblur'}
        % Should be in meters.  The value is set for shift invariant blur.
        % The value for birefringent could come from here, too.
        oi.diffuser.blur = val;

        % Precomputed shift-variant (sv) psf and related parameters          
    case {'psfstruct','shiftvariantstructure'}
        % This structure
        oi.psf = val;
    case {'svpsf','sampledrtpsf','shiftvariantpsf'}
        % The precomputed shift-variant psfs
        oi.psf.psf = val;
    case {'psfanglestep','psfsampleangles'}
        % Vector of sample angles  
        oi.psf.sampAngles= val;
    case {'psfopticsname','raytraceopticsname'}
        % Name of the optics data are derived from
        oi.psf.opticsName =val;
    case 'psfimageheights'
        % Vector of sample image heights
        oi.psf.imgHeight = val;
    case 'psfwavelength'
        % Wavelengths for this calculation. Should match the optics, I
        % think.  Not sure why it is duplicated.
        oi.psf.wavelength = val;

    case 'depthmap'
        % Depth map, usuaully inherited from scene, in meters
        % oiSet(oi,'depth map',dMap);
        oi.depthMap = val;
        
    otherwise
        error('Unknown oiSet parameter: %s',parm);
end

end