function scene = sceneSet(scene,parm,val,varargin)
%Set ISET scene parameter values
%
%   scene = sceneSet(scene,parm,val,varargin)
%
%  All of the parameters of a scene structure are set through the calls to
%  this routine.
%
%  The scene is the object; parm is the name of the parameter; val is the
%  value of the parameter; varargin allows some additional parameters in
%  certain cases.
%
%  There is a corresponding sceneGet routine.  Many fewer parameters are
%  available for 'sceneSet' than 'sceneGet'. This is because many of the
%  parameters derived from sceneGet are derived from the few parameters
%  that can be set, and sometimes the derived quantities require some
%  knowledge of the optics as well.
%
% Parameters are freely formatted with spaces and capitals, which are
% removed in this routine. 
%
%  Examples:
%    scene = sceneSet(scene,'name','myScene');   % Set the scene name
%    scene = sceneSet(scene,'name','myScene');   % Set the scene name
%    scene = sceneSet(scene,'h fov',3);          % Set scene horizontal field of view to 3 deg
% 
%
% Scene parameters - 
%
%  Basic
%      {'name'}          - An informative name describing the scene
%      {'type'}          - The string 'scene'
%      {'distance'}      - Object distance from the optics (meters)
%      {'wangular'}      - Width (horizontal) field of view
%      {'magnification'} - Always 1 for scenes.
%
% Scene radiance parameters
%      {data}   - structure containing the data
%        {'photons' }
%         row x col x nWave array representing the scene radiance in
%         photons.  Typically these values are only stored with single
%         precision to save space.
%         If you set bit depth to 64, they will be stored as double
%               s = sceneSet(s,'bit depth',64);
%
%        {'peak photon radiance'} - Used for monochromatic scenes mainly; not a
%            variable, but a function
%
%         N.B. After writing to the 'photons' field we clear the luminance
%         data.  To fill the luminance slot use:
%
%            [lum, meanL] = sceneCalculateLuminance(scene);
%            scene = sceneSet(scene,'luminance',lum);
%
%         N.B. To adjust the scene mean luminance, use sceneAdjustLuminance
% Depth
%      {'depth map'} - Stored in meters.  Used with RenderToolbox
%      synthetic scenes.  (See scene3D pdcproject directory).
%
% Scene color information
%      {'spectrum'}   - structure that contains wavelength information
%        {'wavelength'} - Wavelength sample values (nm)
%
% Illuminant
%     {'illuminant'}  - Scene illumination structure
%      {'illuminant energy'}  - Illuminant spd in energy is stored W/sr/nm/sec
%      {'illuminant photons'} - Photons are converted to energy and stored 
%      {'illuminant comment'} - Comment
%      {'illuminant name'}    - Identifier for illuminant.
%         See sceneIlluminantScale() for setting the illuminant level in
%         certain cases of unknown reflectance and illuminant conditions.
%         Though, this is being deprecated.
%
% Private variables used by ISET but not set by the user
%
%    Used for management of compressed photons
%      {'bitdepth'}
%      {'knownReflectance'} - For scenes when a reflectance is known
%                             (reflectance,i,j,w)
%
%    Used to cache the scene luminance rather than recompute
%      {'luminance'}
%      {'meanluminance'}
%
%    Auxiliary
%      {'consistency'} - Display consistent with window data
%
% Copyright ImagEval Consultants, LLC, 2003.

if ~exist('parm','var') || isempty(parm), error('Param must be defined.'); end
if ~exist('val','var'), error('Value field required.'); end  % empty is OK

parm = ieParamFormat(parm);

switch parm 
    case {'name','scenename'}
        if ischar(val), scene.name = val; else error('val is not str'); end
    case 'type'
        scene.type = val;
    case {'filename'}
        % When the data are ready from a file, we may save the file name.
        % Happens, perhaps, when reading multispectral image data.
        % Infrequently used.
        scene.filename = val;
    case {'consistency','computationalconsistency'}
        % When parameters are changed, the consistency flag on the optical
        % image changes.  This is irrelevant for the scene case.
        scene.consistency = val;
    case {'gamma'}
        % sceneSet([],'gamma',1);
        % Should this be ieSessionSet('scene gamma',val)
        % hObj = ieSessionGet('scene window');
        hdl = ieSessionGet('scene window handle');
        % eventdata = [];
        set(hdl.editGamma,'string',num2str(val));
        % sceneWindow('sceneRefresh',hObj,eventdata,hdl);
        sceneWindow;
        
    case {'distance' }
        % Positive for scenes, negative for optical images
        scene.distance = val;

    case {'wangular','widthangular','hfov','horizontalfieldofview','fov'}
        if val > 180, val = 180 - eps; warndlg('Warning: fov > 180');
        elseif val < 0, val = eps; warndlg('Warning fov < 0');
        end
        scene.wAngular = val;

    case 'magnification'
        % Scenes should always have a magnification of 1.
        if val ~= 1, warndlg('Scene must have magnification 1'); end
        scene.magnification = 1;

    case {'data','datastructure'}
        scene.data = val;
        
    case {'photons','cphotons','compressedphotons'}
        % sceneSet(scene,'photons',val,[wavelength])
        % val is typically a 3D (row,col,wave) matrix.
      
        bitDepth = sceneGet(scene, 'bitDepth');
        if isempty(bitDepth), error('Compression parameters required'); end
        if ~isempty(varargin)
            idx = ieFindWaveIndex(sceneGet(scene, 'wave'), varargin{1});
            idx = logical(idx);
        end
        
        switch bitDepth
            case 64 % Double
                if isempty(varargin)
                    scene.data.photons = val;
                else
                    scene.data.photons(:,:,idx) = val;
                end
            case 32 % Single
                if isempty(varargin)
                    scene.data.photons = single(val);
                else
                    scene.data.photons(:,:,idx) = single(val);
                end
            otherwise
                error('Unsupported data bit depth %f\n',bitDepth);
        end

        % Clear out derivative luminance/illuminance computations
        scene = sceneSet(scene, 'luminance', []);

    case 'energy'
        % scene = sceneSet(scene,'energy',energy,wave);
        % 
        % The user specified the scene in units of energy.  We convert to
        % photons and set the data as compressed photons for them.
        %
        wave = sceneGet(scene,'wave');
        photons = zeros(size(val));
        [r,c,w] = size(photons);
        if w ~= length(wave), error('Data mismatch'); end
        
        % h = waitbar(0,'Energy to photons');
        for ii=1:w
            % waitbar(ii/w,h);           
            % Get the first image plane from the energy hypercube.
            % Make it a row vector
            tmp = val(:,:,ii); tmp = tmp(:)';
            % Convert the rwo vector from energy to photons
            tmp = Energy2Quanta(wave(w),tmp);
            % Reshape it and place it in the photon hypercube
            photons(:,:,ii) = reshape(tmp,r,c);
        end
        % close(h);
        scene = sceneSet(scene,'cphotons',photons);

    case {'peakradiance','peakphotonradiance'}
        % Used with monochromatic scenes to set the radiance in photons.
        % scene = sceneSet(scene,'peak radiance',1e17);
        oldPeak = sceneGet(scene,'peak radiance');
        p  = sceneGet(scene,'photons');
        scene = sceneSet(scene,'cphotons',val*(p/oldPeak));
    case {'peakenergyradiance'}
        % As above, but for energy.  Useful for equating energy in a series
        % of monochromatic images.
        error('Peak energy radiance not yet implemented');
    case {'depthmap'}
        % Depth map is always in meters
        scene.depthMap = val;
        
        % As of July 2014 or so, these are no longer used
    case {'datamin','dmin'}
        % These are photons (radiance)
        scene.data.dmin = val;
    case {'datamax','dmax'}
        % These are photon (radiance)
        scene.data.dmax = val;
        %%%%%%%%%%%%%%
        
    case 'bitdepth'
        % Bit depth controls whether the data are stored as single (32) or
        % double (64)
        if val ~= 32 && val ~=64, error('Bad bit depth %i\n',val);
        else                      scene.data.bitDepth = val;
        end
        % scene = sceneClearData(scene);
        
        % Not sure this is used much or at all any more.  It is in
        % sceneIlluminantScale alone, as far as I can tell. - BW
    case 'knownreflectance'
        % We  store a known reflectance at location (i,j) for wavelength
        % w. This information is used to set the illuminant level properly
        % and to keep track of reflectances.
        if length(val) ~= 4 || val(1) > 1 || val(1) < 0
            error('known reflectance is [reflectance,row,col,wave]'); 
        end
        scene.data.knownReflectance = val;

    case {'luminance','lum'}
        % The value here is stored to make computation efficient.  But it
        % is dangerous because this value could be inconsistent with the
        % photons if we are not careful.
        if strcmp(sceneGet(scene,'type'),'scene'), scene.data.luminance = val;
        else error('Cannot set luminance of a non-scene structure.');
        end
    case {'meanluminance','meanl'}
        % This leaves open the possibility that the mean differs from the
        % mean calculated from the real luminance data.  We should probably
        % have this set by a sceneAdjustLuminance() call.
        scene = sceneAdjustLuminance(scene,val);
        scene.data.meanL = val;
        % Get this working
    case {'spectrum','wavespectrum','wavelengthspectrumstructure'}
        scene.spectrum  = val;
        %     case {'binwidth','wavelengthspacing'}
        %         scene.spectrum.binwidth = val;
    case {'wave','wavelength','wavelengthnanometers'}
        % scene = sceneSet(scene,'wave',wave);
        % If there are data, we interpolate the data as well as setting the
        % wavelength.
        % If there are no data, we just set the wavelength.

        if ~checkfields(scene,'data','photons') || isempty(scene.data.photons)
            % No data, so just set the spectrum
            scene.spectrum.wave = val;
        else
            % Because there are data present, we must interpolate the
            % photon data
            scene = sceneInterpolateW(scene,val);
        end
        
        % Scene illumination information
    case {'illuminant'}
        % The whole structure
        scene.illuminant = val;
    case {'illuminantdata','illuminantenergy'}
        % This set changes the illuminant, but it does not change the
        % radiance SPD.  Hence, changing the illuminant (implicitly)
        % changes the reflectance. This might not be what you want.  If you
        % want to change the scene as if it is illuminanted differently,
        % use the function: sceneAdjustIlluminant()
        
        % The data are stored in energy  units, unfortunately The data can
        % be a vector (one SPD for the whole image) or they can be a RGB
        % format SPD with a different illuminant at each position.
        illuminant = sceneGet(scene,'illuminant');
        illuminant = illuminantSet(illuminant,'energy',val);
        scene = sceneSet(scene,'illuminant',illuminant);
        % scene.illuminant = illuminantSet(scene.illuminant,'energy',val(:));
    case {'illuminantphotons'}
        % See comment above about sceneAdjustIlluminant.
        %
        % sceneSet(scene,'illuminant photons',data)
        %
        % We have to handle the spectral and the spatial spectral cases
        % within the illuminantSet.  At this point, val can be a vector or
        % an RGB format matrix.
        if checkfields(scene,'illuminant')
            scene.illuminant = illuminantSet(scene.illuminant,'photons',val);
        else
            % We use a default d65.  The user must change to be consistent
            wave = sceneGet(scene,'wave');
            scene.illuminant = illuminantCreate('d65',wave);
            scene.illuminant = illuminantSet(scene.illuminant,'photons',val);
        end
    case {'illuminantname'}
        scene.illuminant = illuminantSet(scene.illuminant,'name',val);
    case {'illuminantwave'}
        error('Call scene set wave, not illuminant wave');
    case {'illuminantcomment'}
        scene.illuminant.comment = val;
    case {'illuminantspectrum'}
        scene.illuminant.spectrum = val;
        
    otherwise
        disp(['Unknown sceneSet parameter: ',parm]);
end

end
