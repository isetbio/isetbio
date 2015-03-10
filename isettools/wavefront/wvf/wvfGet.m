function val = wvfGet(wvf,parm,varargin)
% Get wavefront structure parameters and derived properties
%
%    val = wvfGet(wvf,parm,varargin)
%
% Wavefront properties are either stored as parameters or computed from those
% parameters. We generally store only unique values and  calculate all
% derived values.
%
%  A '*' indicates that the syntax wvfGet(wvf,param,unit) can be used, where
%  unit specifies the spatial scale of the returned value:
%    length: 'm', 'cm', 'mm','um', 'nm'.
%    angle: 'deg', 'min', 'sec'
%
%  A leading '+ indicates that this is a get only parameter and may not be
%  set.
%
% Parameters:
%
%  Bookkeeping
%   'name' - Name of this object
%   'type' - Type of this object, should always be 'wvf'
%
%  Measured parameter and Zernike coefficients 
%   'zcoeffs' - Zernike coefficients, , OSA standard numbering/coords
%   'measured pupil size' - Pupil size for wavefront aberration meaurements (mm,*)
%   'measured wl' - Wavefront aberration measurement wavelength (nm,*)
%   'measured optical axis' - Measured optical axis (deg)
%   'measured observer accommodation' - Observer accommodation at aberration measurement time (diopters)
%   'measured observer focus correction' - Focus correction added optically for observer at measurement time (diopters)
%
%  Spatial sampling parameters
%    'sample interval domain' - Which domain has sample interval held constant with wavelength ('psf', 'pupil')
%    'number spatial samples' - Number of spatial samples (pixel) for pupil function and psf
%    'ref pupil plane size' - Size of sampled pupil plane at measurement wavelength (mm,*)
%    'ref pupil plane sample interval' - Pixel sample interval in pupil plane at measurement wavelength (mm,*)
%    'ref psf sample interval' - Sampling interval for psf at measurment wavelength (arcminute/pixel)
%  + 'pupil plane size' - Size of sampled pupil plane at any calculated wavelength(s) (mm)
%  + 'psf arcmin per sample' - Sampling interval for psf at any calculated wavelength(s) (min/pixel)
%  + 'psf angle per sample' - Sampling interval for psf at any calculated wavelength(s) (min,*/pixel)
%  + 'psf angular samples' - One-d slice of sampled angles for psf, centered on 0, for a single wavelength (min,*)
%  + 'psf spatial samples' - One-d slice of sampled psf in spatial units, centered on 0 for a single wavelength (*)
%  + 'pupil spatial samples' - One-d slice of sampled pupil function in spatial units, centered on 0 for a single wavelength (*)
%  + 'middle row' - The middle row of sampled functions
%
%  Calculation parameters
%     'calc pupil size'  - Pupil size for calculation (mm,*)
%     'calc optical axis' - Optical axis to compute for (deg)
%     'calc observer accommodation' - Observer accommodation at calculation time (diopters)
%     'calc observer focus correction' - Focus correction added optically for observer at calculation time (diopters)
%     'calc wavelengths' - Wavelengths to calculate over (nm,*)
%     'calc cone psf info' - Structure with cone sensitivities and weighting spectrum for computing cone psfs.
%  +  'number calc wavelengths' - Number of wavelengths to calculate over
%
% Pupil and sointspread function
%  +  'wavefront aberrations' - The wavefront aberrations in microns.  Must call wvfComputePupilFunction on wvf before get (um)
%  +  'pupil function' - The pupil function.  Must call wvfComputePupilFunction on wvf before get.
%  +  'psf' - Point spread function.  Must call wvfComputePSF on wvf before get
%  +  'psf centered' - Peak of PSF is at center of returned matrix
%  +  '1d psf' - One dimensional horizontal (along row) slice through PSF centered on its max
%  +  'diffraction psf' - Diffraction limite PSF
%  +  'cone psf' - PSF as seen by cones for given weighting spectrum.
%
% Stiles Crawford Effect
%     'sce params' - The whole structure
%     'sce x0'
%     'sce y0'
%     'sce rho'
%     'sce wavelengths'*
%  +  'sce fraction' - How much light is effectively lost by cones because of sce
%  +  'areapix' - Used in computation of sce fraction
%  +  'areapixapod' - Used in computation of sce fraction
%  +  'cone sce fraction' - SCE fraction for cone psfs
%
% Need to be implemented/checked/documented
%  +  'distanceperpix'
%  +  'samplesspace'
%  +  'strehl'     - Ratio of peak of diffraction limited to actual
%
% Examples:
% * Compute diffraction limited psf
%   wvfP = wvfCreate;
%   wvfP = wvfComputePSF(wvfP);
%   vcNewGraphWin; wvfPlot(wvfP,'image psf','um',550);
%
%   psf = wvfGet(wvfP,'diffraction psf',550); vcNewGraphWin; mesh(psf)
%
% * Strehl is ratio of diffraction and current
%   wvfP = wvfComputePSF(wvfP); wvfGet(wvfP,'strehl',550)
%
% * Blur and recompute.  4th coefficient is defocus
%   z = wvfGet(wvfP,'zcoeffs');z(4) = 0.3; wvfP = wvfSet(wvfP,'zcoeffs',z);
%   wvfP = wvfComputePSF(wvfP); wvfGet(wvfP,'strehl',550)
%
%   wvf = wvfCreate; wvf = wvfComputePSF(wvf); 
%   otf = wvfGet(wvf,'otf'); f = wvfGet(wvf,'otf support','mm');
%   vcNewGraphWin; mesh(f,f,otf);
%
%   lsf = wvfGet(wvf,'lsf'); x = wvfGet(wvf,'lsf support','mm');
%   vcNewGraphWin; plot(x,lsf);
%
% See also: wvfSet, wvfCreate, wvfComputePupilFunction, wvfComputePSF,
% sceCreate, sceGet
%
% (c) Wavefront Toolbox Team 2011, 2012

if ~exist('parm','var') || isempty(parm), error('Parameter must be defined.'); end

% Default is empty when the parameter is not yet defined.
val = [];

parm = ieParamFormat(parm);
if strcmp(parm,'wave'), warning('Change wave to calc wave'); dbstop; end
if strcmp(parm,'samplesamples'), warning('Change wave to nsamples'); dbstop; end

%% We will subdivide the gets over time
%  We plan to create get functions, such as wvfpsfGet(), or wvfsceGet, to
%  introduce some more order.  We will modify ieParameterOtype to help with
%  this.
%
switch parm
    %% Book-keeping
    case 'name'
        val = wvf.name;
    case 'type'
        val = wvf.type; 
    
        %% Pupil plane properties
        %
        % The Zernicke coefficients define the wavefront aberrations in the
        % pupil plane.  Various quantities are derived from this.
        %
        % This group contains many parameters related to the pupil
        % functions
    case {'zcoeffs','zcoeff','zcoef'}
        % Zernike coeffs
        % wvfGet(wvf,'zcoeffs',idx);
        % idx is optional, and can be a vector of j values
        % or a string array of coefficient names (see wvfOSAIndexToVectorIndex).
        % Note that j values start at 0, and that is the convention followed
        % here.  If idx is passed, the length of val matches that of idx.
        % And, it is an error if you try to get a coefficient that has not
        % been set.
        if (isempty(varargin))
            val = wvf.zcoeffs;
        else
            idx = wvfOSAIndexToVectorIndex(varargin{1});
            tempcoeffs = wvf.zcoeffs;
            maxidx = max(idx);
            if (maxidx > length(wvf.zcoeffs))
                tempcoeffs(length(tempcoeffs)+1:maxidx) = 0;
            end
            val = tempcoeffs(idx);
        end
       case {'wavefrontaberrations'}
        % The wavefront aberrations are derived from Zernicke coefficients
        % in the routine wvfComputePupilFunction
        % 
        % If there are multiple wavelengths, then this is a cell array of
        % matrices wvfGet(wvf,'wavefront aberrations',wList) This comes
        % back in microns, and if I were a better person I would have
        % provided unit passing and conversion.
        
        % Can't do the get unless it has already been computed and is not stale.
        if (~isfield(wvf,'pupilfunc') || ...
                ~isfield(wvf,'PUPILFUNCTION_STALE') || ...
                wvf.PUPILFUNCTION_STALE)
            error('Must compute wavefront aberrations before getting them.  Use wvfComputePupilFunction or wvfComputePSF.');
        end
        
        % Return whole cell array of wavefront aberrations over wavelength if
        % no argument passed.  If there is just one wavelength, we
        % return the .wavefront aberrations as a matrix, rather than as a cell
        % array with one entry.
        if isempty(varargin)
            if (length(wvf.wavefrontaberrations) == 1)
                val = wvf.wavefrontaberrations{1};
            else
                val = wvf.wavefrontaberrations;
            end
        else
            wList = varargin{1}; idx = wvfWave2idx(wvf,wList);
            nWave = wvfGet(wvf,'nwave');
            if idx > nWave, error('idx (%d) > nWave',idx,nWave);
            else val = wvf.wavefrontaberrations{idx};
            end
        end
        
    case {'pupilfunction','pupilfunc','pupfun'}
        % The pupil function is derived from Zernicke coefficients in the
        % routine wvfComputePupilFunction 
        %
        % If there are multiple wavelengths, then this is a cell array of
        % matrices
        %   wvfGet(wvf,'pupilfunc',wList)
        
        % Can't do the get unless it has already been computed and is not
        % stale.
        if (~isfield(wvf,'pupilfunc') || ...
                ~isfield(wvf,'PUPILFUNCTION_STALE') || ...
                wvf.PUPILFUNCTION_STALE)
            error('Must compute pupil function before getting it.  Use wvfComputePupilFunction or wvfComputePSF.');
        end
        
        % Return whole cell array of pupil functions over wavelength if
        % no argument passed.  If there is just one wavelength, we
        % return the pupil function as a matrix, rather than as a cell
        % array with one entry.
        if isempty(varargin)
            if (length(wvf.pupilfunc) == 1)
                val = wvf.pupilfunc{1};
            else
                val = wvf.pupilfunc;
            end
        else
            wList = varargin{1}; idx = wvfWave2idx(wvf,wList);
            nWave = wvfGet(wvf,'nwave');
            if idx > nWave, error('idx (%d) > nWave',idx,nWave);
            else val = wvf.pupilfunc{idx};
            end
        end
        
            
        %% The set of measured properties 
        % These form the backdrop for the calculation parameters. 
        %
    case {'measuredpupildiameter','pupilsizemeasured','measuredpupilsize', 'measuredpupil', 'measuredpupilmm'}
        % Pupil diameter in mm over for which wavefront expansion is valid
        % wvfGet(wvf,'measured pupil','mm')
        % wvfGet(wvf,'measured pupil')
        val = wvf.measpupilMM;
        if ~isempty(varargin)
            % Convert to meters and then scale
            val = (val*1e-3)*ieUnitScaleFactor(varargin{1});
        end

    case {'measuredwavelength','wlmeasured','wavelengthmeasured','measuredwl'}
        % Measurement wavelength (nm)
        val = wvf.measWlNM;
        if ~isempty(varargin)
            % Convert to meters and then scale
            val = (val*1e-9)*ieUnitScaleFactor(varargin{1});
        end
        
    case {'measuredopticalaxis','opticalaxismeasued', 'measuredopticalaxisdeg'}
        % Measurement optical axis, degrees eccentric from fovea
        val = wvf.measOpticalAxisDeg;
        
    case {'measuredobserveraccommodation', 'measuredobserveraccommodationdiopters'}
        % Observer accommodation, in diopters relative to relaxed state of eye
        val = wvf.measObserverAcommodationDiopters;
        
    case {'measuredobserverfocuscorrection', 'measuredobserverfocuscorrectiondiopters'}
        % Focus correction added optically for observer at measurement time (diopters)
        val = wvf.measObserverAcommodationDiopters;
        
        %% Spatial sampling parameters related to ...
        %
        % Say more here
    case {'sampleintervaldomain'}
        % What's held constant with calculated wavelength.
        % Choices are 'psf' and 'pupil'
        % This really needs a better explanation.  It has to do with
        % accounting for the index of refraction, undoubtedly.
        val = wvf.constantSampleIntervalDomain;
        
    case {'numberspatialsamples','nsamples','spatialsamples', 'npixels', 'fieldsizepixels'}
        % Number of pixels for both the pupil and psf planes
        % discretization This is a master value - which means that this is
        % the finest resolution.  
        % Why are there both psf and pupil plane spatial samples?
        % Something about the index of refraction for the separation, but
        % not for the number ...
        val = wvf.nSpatialSamples;
        
    case {'refpupilplanesize', 'refpupilplanesizemm', 'fieldsizemm'}
        % Total size of computed field in pupil plane.  This is for the measurement
        % wavelength and sets the scale for calculations at other
        % wavelengths.  
        %Shouldn't this have 'measured' in the title?
        val = wvf.refSizeOfFieldMM;
        if ~isempty(varargin)
            val = (val*1e-3)*ieUnitScaleFactor(varargin{1});
        end
        
    case {'refpupilplanesampleinterval', 'refpupilplanesampleintervalmm', 'fieldsamplesize', 'fieldsamplesizemmperpixel'}
        % Pixel sample interval of sample pupil field. This is for the measurement
        % wavelength and sets the scale for calculations at other
        % wavelengths.  
        % Shouldn't this have measured in the title?
        val = wvf.refSizeOfFieldMM/wvf.nSpatialSamples;
        if ~isempty(varargin)
            val = (val*1e-3)*ieUnitScaleFactor(varargin{1});
        end
        
    case {'refpsfsampleinterval' 'refpsfarcminpersample', 'refpsfarcminperpixel'}
        % Arc minutes per pixel of the sampled psf at the measurement
        % wavelength.  This is for the measurement
        % wavelength and sets the scale for calculations at other
        % wavelengths.
        radiansPerPixel = wvfGet(wvf,'measured wl','mm')/wvfGet(wvf,'ref pupil plane size','mm');
        val = (180*60/3.1416)*radiansPerPixel;
        
         
        %% Calculation parameters
        % The calculation can take place at different wavelengths and pupil
        % diameters than the measurement.  The settings for the calculation
        % are below here, I think.  These should have calc in the title, I
        % think.
        %
        case {'pupilplanesize', 'pupilplanesizemm'}
        % wvfGet(wvf,'pupil plane size',units,wList)
        % Total size of computed field in pupil plane, for calculated
        % wavelengths(s)
        
        % Get wavelengths.  What if varargin{2} is empty?
        wList = varargin{2};
        waveIdx = wvfWave2idx(wvf,wList);
        wavelengths = wvfGet(wvf,'calc wavelengths','nm');
        
        % Figure out what's being held constant with wavelength and act
        % appropriately.
        whichDomain = wvfGet(wvf,'sample interval domain');
        if (strcmp(whichDomain,'psf'))
            val = wvfGet(wvf,'ref pupil plane size','mm')*wavelengths(waveIdx)/wvfGet(wvf,'measured wl','nm');
        elseif (strcmp(whichDomain,'pupil'))
            val = wvfGet(wvf,'ref pupil plane size','mm')*ones(length(waveIdx),1);
        else
            error('Unknown sample interval domain ''%s''',whichDomain);
        end
        
        % Unit conversion.  If varargin{1} is empty, then the units are
        % 'mm' and we leave it alone.
        if ~isempty(varargin)
            val = (val*1e-3)*ieUnitScaleFactor(varargin{1});
        end
        
    case {'calcpupildiameter','calcpupilsize', 'calculatedpupil'}
        % Pupil diameter to use when computing pupil function and PSF.  
        % The calc pupil diameter must
        % be less than or equal to measured pupil size.
        %  wvfGet(wvf,'calculated pupil','mm')
        %  wvfGet(wvf,'calculated pupil','um')
        val = wvf.calcpupilMM;
        
        % Adjust units
        if ~isempty(varargin)
            val = (val*1e-3)*ieUnitScaleFactor(varargin{1});
        end
        
    case {'calcopticalaxis'}
        % Specify optical axis at calculation time
        val = wvf.calcOpticalAxisDegrees;
        if (val ~= wvfGet(wvf,'measuredobserveraccommodation'))
            error('We do not currently know how to deal with values that differ from measurement time');
        end
        
    case {'calcobserveraccommodation'}
        % Specify observer accommodation at calculation time
        val = wvf.calcObserverAccommodationDiopters;
        if (val ~= wvfGet(wvf,'measuredobserveraccommodation'))
            error('We do not currently know how to deal with values that differ from measurement time');
        end
        
    case {'calcobserverfocuscorrection', 'defocusdiopters'}
        % Specify optical correction added to observer focus at calculation time
        val = wvf.calcObserverFocusCorrectionDiopters;
        
    case {'calcwave','calcwavelengths','wavelengths','wavelength','wls'}
        % Wavelengths to compute on
        % wvfGet(wvf,'wave',unit,idx)
        % wvfGet(wvf,'wave','um',3)
        % May be a vector or single wavelength
        val = wvf.wls;
        
        % Adjust units
        if ~isempty(varargin)
            unit = varargin{1};
            val = val*(1e-9)*ieUnitScaleFactor(unit);
        end
        
        % Select wavelength if indices were passed
        if length(varargin) > 1, val = val(varargin{2}); end
        
    case {'calcconepsfinfo'}
        % Weighting spectrum used in calculation of polychromatic psf
        val = wvf.conePsfInfo;
        
    case {'calcnwave','nwave','numbercalcwavelengths','nwavelengths'}
        % Number of wavelengths to calculate at
        val = length(wvf.wls);
       
        %% Point spread parameters
        %  The point spread is an important calculation.
        %  We need linespread and otf, too.
        %
        case 'psf'
        % Get the PSF.
        %   wvfGet(wvf,'psf',wList)
        
        % Force user to code to explicitly compute the PSF if it isn't done.  Not ideal
        % but should be OK.
        if (~isfield(wvf,'psf') || ~isfield(wvf,'PSF_STALE') || wvf.PSF_STALE)
            error('Must explicitly compute PSF on wvf structure before getting it.  Use wvfComputePSF');
        end
        
        % Return whole cell array of psfs over wavelength if
        % no argument passed.  If there is just one wavelength, we
        % return the pupil function as a matrix, rather than as a cell
        % array with one entry.
        if isempty(varargin)
            if (length(wvf.psf) == 1)
                val = wvf.psf{1};
            else
                val = wvf.psf;
            end
        else
            wList = varargin{1}; idx = wvfWave2idx(wvf,wList);
            nWave = wvfGet(wvf,'nwave');
            if idx > nWave, error('idx (%d) > nWave',idx,nWave);
            else val = wvf.psf{idx};
            end
        end
        
    case 'diffractionpsf'
        % Compute and return diffraction limited PSF for the calculation
        % wavelength  and pupil diameter.
        %
        %   wvfGet(wvf,'diffraction psf',wList);
        if ~isempty(varargin), wList= varargin{1};
        else                   wList = wvfGet(wvf,'calc wave');
        end
        zcoeffs = 0;
        wvfTemp = wvfSet(wvf,'zcoeffs',zcoeffs);
        wvfTemp = wvfSet(wvfTemp,'wave',wList(1));
        wvfTemp = wvfComputePSF(wvfTemp);
        val = wvfGet(wvfTemp,'psf',wList(1));
        
    case {'psfarcminpersample', 'psfarcminperpixel', 'arcminperpix'}
        % wvfGet(wvf,'psf arcmin per sample',wList)
        %
        % Arc minutes per pixel in the psf domain, for the calculated
        % wavelength(s).
        
        % Get wavelengths
        wavelengths = wvfGet(wvf,'calc wavelengths','mm');
        wList = varargin{1};
        waveIdx = wvfWave2idx(wvf,wList);
        
        % Figure out what's being held constant with wavelength and act
        % appropriately.
        whichDomain = wvfGet(wvf,'sample interval domain');
        if (strcmp(whichDomain,'psf'))
            val = wvfGet(wvf,'ref psf arcmin per pixel')*ones(length(waveIdx),1);
        elseif (strcmp(whichDomain,'pupil'))
            radiansPerPixel = ...
                wavelengths(waveIdx)/wvfGet(wvf,'ref pupil plane size','mm');
            val = (180*60/pi)*radiansPerPixel;
        else
            error('Unknown sample interval domain ''%s''',whichDomain);
        end
        
    case {'psfanglepersample','angleperpixel','angperpix'}
        % Angular extent per pixel in the psf domain, for calculated
        % wavelength(s).
        %
        % wvfGet(wvf,'psf angle per sample',unit,wList)
        % unit = 'min' (default), 'deg', or 'sec'
        unit  = varargin{1};
        wList = varargin{2};
        val = wvfGet(wvf,'psf arcmin per sample',wList);
        if ~isempty(unit)
            unit = lower(unit);
            switch unit
                case 'deg'
                    val = val/60;
                case 'sec'
                    val = val*60;
                case 'min'
                    % Default
                otherwise
                    error('unknown angle unit %s\n',unit);
            end
        end
        
    case {'psfangularsamples'} % 'samplesangle' 'samplesarcmin','supportarcmin'
        % Return one-d slice of sampled angles for psf, centered on 0, for
        % a single wavelength
        % wvfGet(wvf,'psf angular samples',unit,waveIdx)
        % unit = 'min' (default), 'deg', or 'sec'
        % Should call routine below to get anglePerPix.
        unit = 'min'; wList = wvfGet(wvf,'measured wavelength');
        if ~isempty(varargin), unit = varargin{1}; end
        if (length(varargin) > 1), wList = varargin{2}; end
        if length(wList) > 1
            error('This only works for one wavelength at a time');
        end
        
        anglePerPix = wvfGet(wvf,'psf angle per sample',unit,wList);
        middleRow = wvfGet(wvf,'middle row');
        nPixels = wvfGet(wvf,'spatial samples');
        val = anglePerPix*((1:nPixels)-middleRow);
        
    case {'psfangularsample'}
        % wvfGet(wvf,'psf angular sample',unit,waveIdx)
        % unit = 'min' (default), 'deg', or 'sec'
        unit = varargin{1};
        wList = varargin{2};
        if (length(wList) > 1)
            error('This only works for one wavelength at a time');
        end
        val = wvfGet(wvf,'psf angle per sample',unit,wList);

    case {'psfspatialsamples','samplesspace','supportspace','spatialsupport'}
        % wvfGet(wvf,'samples space','um',wList)
        % Spatial support in samples, centered on 0
        % Unit and wavelength must be specified
        % Should call case below to get one val, and then scale up by
        % number of pixels.
        
        % This parameter matters for the OTF and PSF quite a bit.  It
        % is the number of um per degree on the retina.
        mPerDeg = (330*10^-6);      % Meters per deg
        unit = 'deg'; wave = wvfGet(wvf,'calc wave');
        if ~isempty(varargin), unit = varargin{1}; end
        if length(varargin) > 1, wave = varargin{2}; end
        
        % Get the angular samples in degrees
        val = wvfGet(wvf,'psf angular samples','deg',wave);
        
        % Convert to meters and then to selected spatial scale
        val = val*mPerDeg;  % Sample in meters 
        val = val*ieUnitScaleFactor(unit);
        
    case {'psfspatialsample'}
        % This parameter matters for the OTF and PSF quite a bit.  It
        % is the number of um per degree on the retina.
        umPerDeg = (330*10^-6);
        unit = 'mm'; wList = wvfGet(wvf,'measured wavelength');
        if ~isempty(varargin), unit = varargin{1}; end
        if length(varargin) > 1, wList = varargin{2}; end
        if length(wList) > 1, error('One wavelength only'); end
        
        % Get the samples in degrees
        val = wvfGet(wvf,'psf angular sample','deg',wList);
        
        % Convert to meters and then to selected spatial scale
        val = val*umPerDeg;  % Sample in meters assuming 300 um / deg
        val = val*ieUnitScaleFactor(unit);
        
    case {'pupilspatialsamples'}
        % wvfGet(wvf,'pupil spatial samples','mm',wList)
        % Spatial support in samples, centered on 0
        % Unit and wavelength must be specified
        
        unit = varargin{1}; wList = varargin{2};
        
        % Get the sampling rate in the pupil plane in space per sample
        spacePerSample = wvfGet(wvf,'pupil plane size',unit,wList)/wvfGet(wvf,'spatial samples');
        nSamples = wvfGet(wvf,'spatial samples');
        middleRow = wvfGet(wvf,'middle row');
        val = spacePerSample*((1:nSamples)-middleRow);
        
    case {'pupilspatialsample'}
        % wvfGet(wvf,'pupil spatial sample','mm',wList)
        % Spatial support in samples, centered on 0

        unit = 'mm'; wList = wvfGet(wvf,'calc wave');
        if ~isempty(varargin), unit = varargin{1}; end
        if length(varargin) > 1, wList = varargin{2}; end
        
        % Get the sampling rate in the pupil plane in space per sample
        val = wvfGet(wvf,'pupil plane size',unit,wList)/wvfGet(wvf,'spatial samples');

    case {'middlerow'}
        val = floor(wvfGet(wvf,'npixels')/2) + 1;
        
    case {'otf'}
        % Return the otf from the psf
        %
        % wvfGet(wvf,'otf',wave)
        % What should the units be?  c/(nPSFSamples*units), I think we
        % should deal with this explicitly here.
        
        wave = wvfGet(wvf,'calc wave');
        if ~isempty(varargin), wave = varargin{1}; end
        psf = wvfGet(wvf,'psf',wave);   % vcNewGraphWin; mesh(psf)
        val = fftshift(psf2otf(psf));   % vcNewGraphWin; mesh(val)
        
    case {'otfsupport'}
        % wvfGet(wvf,'otfsupport',unit,wave)
        unit = 'mm'; wave = wvfGet(wvf,'calc wave');
        if ~isempty(varargin),   unit = varargin{1}; end
        if length(varargin) > 1, wave = varargin{2}; end
        
        %         s = wvfGet(wvf,'psf spatial sample',unit,wave);
        %         n = wvfGet(wvf,'nsamples');   % Should specify psf or pupil, but I thik they are the same
        %         val = (0:(n-1))*(s*n);   % Cycles per unit
        %         val = val - mean(val);
        %
        samp = wvfGet(wvf,'samples space',unit,wave);
        nSamp = length(samp);
        dx = samp(2) - samp(1);
        nyquistF = 1 / (2*dx);   % Line pairs (cycles) per unit space
        val = unitFrequencyList(nSamp)*nyquistF;
        

    case {'lsf'}
        % wave = wvfGet(wvf,'calc wave');
        % lsf = wvfGet(wvf,'lsf',unit,wave); vcNewGraphWin; plot(lsf)
        % For the moment, this only runs if we have precomputed the PSF and
        % we have a matching wavelength in the measured and calc
        % wavelengths.  We need to think this through more.
        wave = wvfGet(wvf,'calc wave');
        if length(varargin) > 1, wave = varargin{1}; end
        
        otf  = wvfGet(wvf,'otf',wave);
        mRow = wvfGet(wvf,'middle row');
        val  = fftshift(abs(fft(otf(mRow,:))));
        
    case {'lsfsupport'}
        % wvfGet(wvf,'lsf support');
        %
        unit = 'mm'; wave = wvfGet(wvf,'calc wave');
        if ~isempty(varargin), unit = varargin{1}; end
        if length(varargin) > 1, wave = varargin{2}; end
        val = wvfGet(wvf,'psf spatial samples',unit,wave);
        
        
        %% Stiles Crawford Effect
        % Account for angle sensitivity of the cone photoreceptors
        %
    case 'sce'
        if isfield(wvf,'sceParams'), val = wvf.sceParams; end
        
    case 'scex0'
        if checkfields(wvf,'sceParams','xo'), val = wvf.sceParams.xo;
        else val = 0;
        end
        
    case 'scey0'
        if checkfields(wvf,'sceParams','yo'), val = wvf.sceParams.yo;
        else val = 0;
        end
        
    case {'scewavelength','scewavelengths','scewave'}
        % This returns the wvf wavelength list if there isn't a sceParams
        % structure.  Might be OK.
        % wvfGet(wvf,'sce wavelengths',unit)
        if checkfields(wvf,'sceParams','wavelengths'), val = wvf.sceParams.wavelengths;
        else val = wvf.wls;
        end
        % Adjust units
        if ~isempty(varargin)
            unit = varargin{1};
            val = val*10e-9*ieUnitScaleFactor(unit);
        end
        
    case 'scerho'
        % Get rho value for a particular wavelength
        %  wvfGet(wvf,'rho',waveList)
        if checkfields(wvf,'sceParams','rho'), val = wvf.sceParams.rho;
        else val = zeros(wvfGet(wvf,'nWave'),1);
        end
        
        % Return rho values for selected wavelengths
        if ~isempty(varargin)
            wave = wvfGet(wvf,'sce wave');  % The waves for rho
            wList = varargin{1};
            index = find(ismember(round(wave),round(wList)));
            if ~isempty(index), val = val(index);
            else error('Passed wavelength not contained in sceParams');
            end
        end
        
    case {'scefraction','scefrac','stilescrawfordeffectfraction'}
        % How much light is effectively lost at each wavelength during
        % cone absorption becauseof the Stiles-Crawford effect.  The most likely
        % use of this is via the scefrac get above.
        %
        % This is computed with the pupil function, and is thus stale
        % if the pupil function is stale.
        if (~isfield(wvf,'pupilfunc') || ...
                ~isfield(wvf,'PUPILFUNCTION_STALE') || ...
                wvf.PUPILFUNCTION_STALE)
            error('Must compute pupil function  before retrieving %s. Use wvfComputePupilFunction or wvfComputePSF', parm);
        end
        
        if isempty(varargin)
            val = wvfGet(wvf,'area pixapod') ./ wvfGet(wvf,'areapix');
        else
            wList = varargin{1};
            val = wvfGet(wvf,'area pixapod',wList) ./ wvfGet(wvf,'areapix',wList);
        end
        
    case {'areapix'}
        % This is the summed amplitude of the pupil function *before*
        % Stiles-Crawford correction over the pixels where the pupil
        % function is defined.  It doesn't have much physical significance,
        % but taking the ratio with areapixapod (just below) tells us
        % how much light is effectively lost at each wavelength during
        % cone absorption becauseof the Stiles-Crawford effect.  The most likely
        % use of this is via the scefrac get above.
        %
        % This is computed with the pupil function, and is thus stale
        % if the pupil function is stale.
        if (~isfield(wvf,'pupilfunc') || ...
                ~isfield(wvf,'PUPILFUNCTION_STALE') || ...
                wvf.PUPILFUNCTION_STALE)
            error('Must compute pupil function  before retrieving %s. Use wvfComputePupilFunction or wvfComputePSF', parm);
        end
        
        if isempty(varargin)
            val = wvf.areapix;
        else
            wList = varargin{1}; idx = wvfWave2idx(wvf,wList);
            nWave = wvfGet(wvf,'nwave');
            if idx > nWave, error('idx (%d) > nWave',idx,nWave);
            else val = wvf.areapix(idx);
            end
        end
        
    case {'areapixapod'}
        % This is the summed amplitude of the pupil function *after*
        % Stiles-Crawford correction over the pixels where the pupil
        % function is defined.  It doesn't have much physical significance,
        % but taking the ratio with areapixapod (just above) tells us
        % how much light is effectively lost at each wavelength during
        % cone absorption becauseof the Stiles-Crawford effect.  The most likely
        % use of this is via the scefrac get above.
        %
        % This is computed with the pupil function, and is thus stale
        % if the pupil function is stale.
        if (~isfield(wvf,'pupilfunc') || ~isfield(wvf,'PUPILFUNCTION_STALE') || ...
                wvf.PUPILFUNCTION_STALE)
            error('Must compute pupil function  before retrieving %s. Use wvfComputePupilFunction or wvfComputePSF', parm);
        end
        
        if isempty(varargin)
            val = wvf.areapixapod;
        else
            wList = varargin{1}; idx = wvfWave2idx(wvf,wList);
            nWave = wvfGet(wvf,'nwave');
            if idx > nWave, error('idx (%d) > nWave',idx,nWave);
            else val = wvf.areapixapod(idx);
            end
        end
        
    case {'sceconesfraction','conescefraction'}
        % SCE fraction for cone psfs
        
        % Can't do this unless psf is computed and not stale
        if (~isfield(wvf,'psf') || ~isfield(wvf,'PSF_STALE') || wvf.PSF_STALE)
            error('Must compute PSF on wvf structure before retrieving %s. Use wvfComputePSF', parm);
        end
        
        [nil,val] = wvfComputeConePSF(wvf);     
        
        
    case 'strehl'
        % Strehl ratio. The strehl is the ratio of the peak of diff limited and the
        % existing PSF at each wavelength.
        %   wvfGet(wvf,'strehl',wList);
        
        % Force user to code to explicitly compute the PSF if it isn't done.  Not ideal
        % but should be OK.
        if (~isfield(wvf,'psf') || ~isfield(wvf,'PSF_STALE') || wvf.PSF_STALE)
            error('Must compute PSF on wvf structure before retrieving %s. Use wvfComputePSF', parm);
        end
        
        % We could write this so that with no arguments we return all of
        % the ratios across wavelengths.  For now, force a request for a
        % wavelength index.
        wList = varargin{1};
        psf = wvfGet(wvf,'psf',wList);
        dpsf = wvfGet(wvf,'diffraction psf',wList);
        val = max(psf(:))/max(dpsf(:));
        
        %         areaPixapod = wvfGet(wvf,'area pixapod',waveIdx);
        %         val = max(psf(:))/areaPixapod^2;
        %         % Old calculation was done in the compute pupil function routine.
        % Now, we do it on the fly in here, for a wavelength
        % strehl(wl) = max(max(psf{wl}))./(areapixapod(wl)^2);
        
    case 'psfcentered'
        % PSF entered so that peak is at middle position in coordinate grid
        %   wvfGet(wvf,'psf centered',wList)
        
        % Force user to code to explicitly compute the PSF if it isn't done.  Not ideal
        % but should be OK.
        if (~isfield(wvf,'psf') || ~isfield(wvf,'PSF_STALE') || wvf.PSF_STALE)
            error('Must compute PSF on wvf structure before retrieving %s. Use wvfComputePSF', parm);
        end
        
        if isempty(varargin), wList = wvfGet(wvf,'wave');
        else wList = varargin{1};
        end
        if length(wList) > 1, error('Only one wavelength permitted');
        else                  val = psfCenter(wvfGet(wvf,'psf',wList));
        end
        
    case '1dpsf'
        % One dimensional slice through the PSF.
        %   wvfGet(wvf,'1d psf',wList,row)
        
        % Force user to code to explicitly compute the PSF if it isn't done.  Not ideal
        % but should be OK.
        if (~isfield(wvf,'psf') || ~isfield(wvf,'PSF_STALE') || wvf.PSF_STALE)
            error('Must compute PSF on wvf structure before retrieving %s. Use wvfComputePSF', parm);
        end
        
        % Defaults
        wList = wvfGet(wvf,'calc wave');
        whichRow = wvfGet(wvf,'middle row');
        
        % Override with varargins
        if ~isempty(varargin),   wList    = varargin{1}; end
        if length(varargin) > 1, whichRow = varargin{2}; end
        
        psf = psfCenter(wvfGet(wvf,'psf',wList));
        val = psf(whichRow,:);
        
    case 'conepsf'
        % PSF as seen by cones for specified weighting spectrum
        
        % Force user to code to explicitly compute the PSF if it isn't done.  Not ideal
        % but should be OK.
        if (~isfield(wvf,'psf') || ~isfield(wvf,'PSF_STALE') || wvf.PSF_STALE)
            error('Must compute PSF on wvf structure before retrieving %s. Use wvfComputePSF', parm);
        end
        
        % Defaults
        val = wvfComputeConePSF(wvf);
    otherwise
        error('Unknown parameter %s\n',parm);
        
end



return

