function val = sensorGet(sensor,param,varargin)
%Get data from ISET image sensor array
%
%     val = sensorGet(sensor,param,varargin)
%
% The (very long) sensor parameter list is described below.  The image
% sensory array is often referred to as sensor, or sensor in the code. The
% sensor array includes a pixel data structure, and this structure has its
% own accessor functions.  The pixel optics depend on microlens structures,
% and there is a separate microlens analysis toolkit.
%
% A '*' indicates a routine that can return different spatial units.
%
% The pixel structure is attached to the sensor.  Ordinarily, you can
% address the pixel structure using
%    pixel = sensorGet(sensor,'pixel');
%    v = pixelGet(pixel,'param');
%
% To simplify the code, you can also use sensorGet to retrieve pixel
% parameters 
%    v = sensorGet(sensor,'pixel param');
%
% Examples:
%    val = sensorGet(sensor,'name')
%    val = sensorGet(sensor,'size');          % row,col
%    val = sensorGet(sensor,'dimension','um');
%    val = sensorGet(sensor,'Electrons',2);   % Second color type
%    val = sensorGet(sensor,'fov')            % degrees
%    val = sensorGet(sensor,'PIXEL')
%    val = sensorGet(sensor,'exposureMethod');% Single, bracketed, cfa
%    val = sensorGet(sensor,'nExposures')     % Number of exposures
%    val = sensorGet(sensor,'filtercolornames')
%    val = sensorGet(sensor,'exposurePlane'); % For bracketing simulation
%    val = sensorGet(sensor,'pixel size','um');
%    val = sensorGet(sensor,'pixel voltage swing');
%
% Basic sensor array parameters
%      {'name'}                 - this sensor name
%      {'type'}                 - always 'sensor'
%      {'row'}                  - sensor rows
%      {'col'}                  - sensor columns
%      {'size'}                 - (rows,cols)
%      {'height'}*              - sensor height (units)
%      {'width'}*               - sensor width  (units)
%      {'dimension'}*           - (height,width)
%      {'spatial support'}*      - position of pixels.
%      {'wspatial resolution'}*  - spatial distance between pixels (width)
%      {'hspatial resolution'}*  - spatial distance between pixels (height)
%
%  Field of view and sampling density
%      {'hfov'}   - horizontal field of view (deg)
%      {'vfov'}   - vertical field of view (deg)
%      {'h deg perpixel'} - horizontal deg per pixel
%      {'v deg perpixel'} - vertical deg per pixel
%      {'h deg perdistance'} - deg per unit horizontal distance *
%      {'v deg perdistance'} - deg per unit vertical distance *
%
% * Units are meters by default but can be set to 'um','mm', etc.
%
%  Sensor optics related
%      {'fov'}                     - sensor horizontal field of view
%      {'chief Ray Angle'}         - chief ray angle in rad at each pixel
%      {'chief Ray Angle Degrees'} - chief ray angle in deg at each pixel
%      {'sensor Etendue'}          - optical efficiency at each pixel
%      {'micro Lens'}              - microlens data structure, accessed
%                                    using mlensGet() and mlensSet (optics
%                                    toolbox only)
%
% Sensor array data
%      {'volts'}          - Sensor output in volts
%      {'photons'}        - Sensor output in number of photons
%      {'photon rate'}    - Photons captured per second
%      {'digital values'} - Sensor output in digital units
%      {'electrons'}      - Sensor output in electrons
%         A single color plane can be returned
%         sensorGet(sensor,'electrons',2);
%      {'dv o rvolts'}    - Return either dv if present, otherwise volts
%      {'roi locs'}       - Stored region of interest (roiLocs)
%      {'roi rect'}       - Rect.  Format is [cmin,rmin,width,height]
%      {'roi volts'}      - Volts inside of stored region of interest
%                           If no ROI stored, ask the user to select.
%      {'roi electrons'}  - Electrons inside of ROI, or user selects
%      {'roi volts mean'} - The mean values in each band
%      {'roi electrons mean'} - As above but electrons
%      {'hline volts'}     - Volts along a horizontal line
%      {'hline electrons'} - horizontal line electrons
%      {'vline volts'}     - vertical line volts
%      {'vline electrons'} - vertical line electrons
%      {'response Ratio'}  - Peak data voltage / largest pixel voltage
%
% Sensor roi
%     {'roi'} - rectangle representing current region of interest
%        Additional ROI processing may be inserted in the next few years.
%
% Sensor array electrical processing properties
%      {'analog Gain'}     - A scale factor that divides the sensor voltage
%                            prior to clipping
%      {'analog Offset'}   - Added to the voltage to stay away from zero,
%                            sometimes used to minimize the effects of dark
%                            noise at the low levels. Formula for offset
%                            and gain: (v + analogOffset)/analogGain)
%
%      {'sensor Dynamic Range'} - dynamic range
%      {'quantization}          -  Quantization structre
%        {'nbits'}              - number of bits in quantization method
%        {'max output'}
%        {'quantization lut'}
%        {'quantization method'}
%
% Sensor color filter array and related color properties
%     {'spectrum'}    - structure about spectral information
%       {'wave'}      - wavelength samples
%       {'binwidth'}  - difference between wavelength samples
%       {'nwave'}     - number of wavelength samples
%     {'color'}
%       {'filter transmissivities'} - Filter transmissivity as a function 
%                                     of wavelength
%       {'infrared filter'}         - Normally the IR, but we sometimes put
%                                     other filters, such as macular
%                                     pigment, in the ir slot.
%
%      {'cfa Name'}     - Best guess at conventional CFA architecture name
%      {'filter Names'} - Cell array of filter names. The first letter of
%                         each filter should indicate the filter color see
%                         sensorColorOrder comments for more information
%      {'nfilters'}    - number of color filters
%      {'filter Color Letters'} - A string with each letter being the first
%                                 letter of a color filter; the letters are
%                                 from the list in sensorColorOrder. The
%                                 pattern field(see below) describes their
%                                 position in array.
%      {'filter Color Letters Cell'} -  As above, but returned in a cell 
%                                       array rather than a string
%      {'filter plotcolors'} - one of rgbcmyk for plotting for this filter
%      {'spectral QE'}       - Product of photodetector QE, IR and color
%                              filters, vignetting or pixel fill factor not
%                              included.
%      {'pattern'}           - Matrix that defines the color filter array
%                              pattern; e.g. [1 2; 2 3] if the spectrra are
%                              RGB and the pattern is a conventional Bayer
%                              [r g; g b]
%
% Noise properties
%      {'dsnu sigma'}           - Dark signal nonuniformity (DSNU) in volts
%      {'prnu sigma'}           - Photoresponse nonuniformity (PRNU) in std
%                                 dev percent
%      {'fpn parameters'}       - (dsnusigma,prnusigma)
%      {'dsnu image'}           - Dark signal non uniformity (DSNU) image
%      {'prnu image',}          - Photo response non uniformity (PRNU)
%      {'column fpn'}           - Column (offset,gain) parameters
%      {'column dsnu'}          - The column offset (Volts)
%      {'column prnu'}          - The column gain (std dev in Volts)
%      {'col offset fpnvector'} - The sensor column offset data
%      {'col gain fpnvector'}   - The sensor column gain data
%      {'noise flag'}           - Governs sensorCompute noise calculations 
%                                 0 no noise at all
%                                 1 shot noise, no electronics noise
%                                 2 shot noise and electronics noise
%      {'reuse noise'}         - Use the stored noise seed
%      {'noise seed'}          - Stored noise seed from last run
%
%  The pixel structure
%      {'pixel'}  - pixel structure is complex; accessed using pixelGet();
%
%  Sensor computation parameters
%      {'auto exposure'}   - Auto-exposure flag (0,1)
%      {'exposure time'}   - Exposure time (sec)
%      {unique exptimes'}  - Unique values from the exposure time list
%      {'exposure plane'}  - Select exposure for display when bracketing
%      {'cds'}             - Correlated double-sampling flag
%      {'pixel vignetting'}- Include pixel optical efficiency in
%             sensorCompute.
%             val = 1 Means vignetting only.
%             val = 2 means microlens included. (Microlens shifting NYI).
%             otherwise, skip both vignetting and microlens.
%      {'sensor compute','sensor compute method'}
%         Swap in a sensorCompute routine.  If this is empty, then the
%         standard vcamera\sensor\mySensorCompute routine will be used.
%      {'nsamples perpixel','npixel samples for computing'}
%         Default is 1.  If parameter is not set, we return the default.
%      {'consistency'}
%         If the consistency field is not present, set it false.  This
%         checks whether the parameters and the displayed image are
%         consistent/updated.
%
% Human sensor special case
%    {'human'} - The structure with all human parameters.  Applies only
%                when the name contains the string 'human' in it
%    {'human lens'}                  - underlying lens structure
%    {'human lens transmittance'}    - lens transmittance
%    {'human lens absorption'}       - lens absorptance
%    {'human macular'}               - macular pigment structure
%    {'human macular densities'}     - macular density
%    {'human macular transmittance'} - macular transmittance
%    {'human macular absorption'}    - macular absorption
%    {'human ocular transmittance'}  - totally transmittance for lens and
%                                      macular pigments
%    {'cone type'}                   - K=1, L=2, M=3 or S=4 (K means none)
%    {'human cone density'}          - densities used to generate mosaic
%    {'rSeed'}                       - seed for generating mosaic
%    {'xy'}                          - xy position of cones in the mosaic
%    {'adaptation gain'}             - cone adaptation gain
%    {'adaptation offset'}           - cone adaptation volts
%    {'adapted volts'}               - adapted volts, will have negatives
%    
%    {'time interval'}          - human eye sampling time interval
%    {'total time'}             - total time of eye movement sequence
%    {'em type'}                - eye movement type vector
%    {'em tremor'}              - eye movement tremor structure
%    {'em drift'}               - eye movement drift structure
%    {'em microsaccade'}        - eye movement microsaccade structure
%
% Sensor motion
%       {'sensor movement'}     - A structure of sensor motion information
%       {'movement positions'}  - Nx2 vector of (x,y) positions in deg
%       {'frames per position'} - N vector of exposures per position
%       {'positions x'}         - 1st column (x) of movement positions
%       {'positions y'}         - 2nd column (y) of movement positions
%
% Miscellaneous - Macbeth color checker (MCC)
%   More chart handling is being introduced.  See chart<TAB>
%
%    {'mcc rect handles'}  - Handles for the rectangle selections in an MCC
%    {'mcc corner points'} - Corner points for the MCC chart
%
%    {'rgb'}               - Display image in sensorImageWindow
%
% See also:  sensorSet
%
% Copyright ImagEval Consultants, LLC, 2005.

if ~exist('param', 'var') || isempty(param)
    error('Param must be defined.');
end

% Default return value.
val = [];

% Is this a pixel call.
[oType,param] = ieParameterOtype(param);
switch oType
    case 'pixel'
        pixel = sensor.pixel;
        if isempty(param), val = pixel;
        elseif   isempty(varargin), val = pixelGet(pixel,param);
        else     val = pixelGet(pixel,param,varargin{1});
        end
        return;
    otherwise
end

% Onward, it is a real sensor call
param = ieParamFormat(param);
switch param
    % Descriptive
    case {'name'}
        if checkfields(sensor,'name'), val = sensor.name; end
    case {'type'}
        if checkfields(sensor,'type'), val = sensor.type; end

        % Size and shape
    case {'rows','row'}
        % There should not be a rows/cols field at all, right, unless the
        % data field is empty?
        if checkfields(sensor,'data','volts')
            val = size(sensor.data.volts,1);
            return;
        elseif checkfields(sensor,'rows'), val = sensor.rows;
        end
    case {'cols','col'}
        % We keep rows/cols field at all, right, unless the
        % data field is empty?
        if checkfields(sensor,'data','volts')
            val = size(sensor.data.volts,2);
            return;
        elseif checkfields(sensor,'cols'), val = sensor.cols;
        end
    case {'size','arrayrowcol'}
        % Note:  dimension is (x,y) but size is (row,col)
        val = [sensorGet(sensor,'rows'),sensorGet(sensor,'cols')];

    case {'height','arrayheight'}
        val = sensorGet(sensor,'rows')*sensorGet(sensor,'deltay');
        if ~isempty(varargin), val = val*ieUnitScaleFactor(varargin{1}); end
    case {'width','arraywidth'}
        val = sensorGet(sensor,'cols')*sensorGet(sensor,'deltax');
        if ~isempty(varargin), val = val*ieUnitScaleFactor(varargin{1}); end

    case {'dimension'}
        val = [sensorGet(sensor,'height'), sensorGet(sensor,'width')];
        if ~isempty(varargin), val = val*ieUnitScaleFactor(varargin{1}); end

        % The resolutions also represent the center-to-center spacing of the pixels.
    case {'wspatialresolution','wres','deltax','widthspatialresolution'}
        PIXEL = sensorGet(sensor,'pixel');
        val = pixelGet(PIXEL,'width') + pixelGet(PIXEL,'widthGap');
        if ~isempty(varargin), val = val*ieUnitScaleFactor(varargin{1}); end

    case {'hspatialresolution','hres','deltay','heightspatialresolultion'}
        PIXEL = sensorGet(sensor,'pixel');
        val = pixelGet(PIXEL,'height') + pixelGet(PIXEL,'heightGap');
        if ~isempty(varargin), val = val*ieUnitScaleFactor(varargin{1}); end

    case {'spatialsupport','xyvaluesinmeters'}
        % ss = sensorGet(sensor,'spatialSupport',units)
        nRows = sensorGet(sensor,'rows');
        nCols = sensorGet(sensor,'cols');
        pSize = pixelGet(sensorGet(sensor,'pixel'),'size');
        val.y = linspace(-nRows*pSize(1)/2 + pSize(1)/2, ...
                          nRows*pSize(1)/2 - pSize(1)/2,nRows);
        val.x = linspace(-nCols*pSize(2)/2 + pSize(2)/2, ...
                          nCols*pSize(2)/2 - pSize(2)/2,nCols);
        if ~isempty(varargin)
            val.y = val.y*ieUnitScaleFactor(varargin{1});
            val.x = val.x*ieUnitScaleFactor(varargin{1});
        end

    case {'chiefrayangle','cra','chiefrayangleradians', ...
          'craradians','craradian','chiefrayangleradian'}
        % Return the chief ray angle for each pixel in radians
        % sensorGet(sensor,'chiefRayAngle',sourceFLMeters)
        support = sensorGet(sensor,'spatialSupport');   %Meters

        % Jst flipped .x and .y positions
        [X,Y] = meshgrid(support.x,support.y);
        if isempty(varargin),
            optics = oiGet(vcGetObject('OI'),'optics');
            sourceFL = opticsGet(optics,'focalLength'); % Meters.
        else
            sourceFL = varargin{1};
        end

        % Chief ray angle of every pixel in radians
        val = atan(sqrt(X.^2 + Y.^2)/sourceFL);

    case {'chiefrayangledegrees', 'cradegrees', ...
          'cradegree','chiefrayangledegree'}
        % sensorGet(sensor,'chiefRayAngleDegrees',sourceFL)
        if isempty(varargin),
            optics = oiGet(vcGetObject('OI'),'optics');
            sourceFL = opticsGet(optics,'focalLength'); % Meters.
        else sourceFL = varargin{1};
        end
        val = rad2deg(sensorGet(sensor,'cra',sourceFL));
    case {'etendue','sensoretendue'}
        % The size of etendue etnrie matches the row/col size of the sensor
        % array. The etendue is computed using the chief ray angle at each
        % pixel and properties of the microlens structure. Routines exist
        % for calculating the optimal placement of the microlens
        % (mlRadiance). We store the bare etendue (no microlens) in the
        % vignetting location.  The improvement due to the microlens array
        % can be calculated by sensor.etendue/sensor.data.vignetting. We
        % need to be careful about clearing these fields and data
        % consistency.
        if checkfields(sensor,'etendue'), val = sensor.etendue; end


    case {'voltage','volts'}
        % sensorGet(sensor,'volts',i) gets the ith sensor data in a vector.
        % sensorGet(sensor,'volts') gets all the sensor data in a plane.
        % This syntax applies to most of the voltage/electron/dv gets
        % below.
        %
        if checkfields(sensor,'data','volts'), val = sensor.data.volts; end
        if ~isempty(varargin)
            val = sensorColorData(val,sensor,varargin{1});
        end
    case {'photonrate'}
        % pRate = sensorGet(sensor,'photon rate');
        % Photons captured per second
        val = sensorGet(sensor,'photons')/sensorGet(sensor,'exposure time');
        
    case{'volts2maxratio','responseratio'}
        v = sensorGet(sensor,'volts');
        pixel = sensorGet(sensor,'pixel');
        sm = pixelGet(pixel,'voltageswing');
        val = max(v(:))/sm;
    case {'analoggain','ag'}
        if checkfields(sensor,'analogGain'), val = sensor.analogGain;
        else val = 1;
        end
    case {'analogoffset','ao'}
        if checkfields(sensor,'analogGain'), val = sensor.analogOffset;
        else   val = 0;
        end
    case {'dv','digitalvalue','digitalvalues'}
        if checkfields(sensor,'data','dv'),val = sensor.data.dv; end
        % Pull out a particular color plane
        if ~isempty(varargin) && ~isempty(val)
            val = sensorColorData(val,sensor,varargin{1});
        end

    case {'electron','electrons','photons'}
        % sensorGet(sensor,'electrons');
        % sensorGet(sensor,'electrons',2);
        % This is also used for human case, where we call the data photons,
        % as in photon absorptions.
        pixel = sensorGet(sensor,'pixel');
        val = sensorGet(sensor,'volts')/pixelGet(pixel,'conversionGain');

        % Pull out a particular color plane
        if ~isempty(varargin)
            val = sensorColorData(val,sensor,varargin{1});
        end
        % Electrons are ints
        val = round(val);

    case {'dvorvolts'}
        val = sensorGet(sensor,'dv');
        if isempty(val), val = sensorGet(sensor,'volts'); end

        % Region of interest for data handling (ROI)
    case {'roi','roilocs'}
        % roiLocs = sensorGet(sensor,'roi');
        % This is the default, which is to return the roi as roi locations,
        % an Nx2 matrix or (r,c) values.
        if checkfields(sensor,'roi')
            % The data can be stored as a rect or as roiLocs.
            val = sensor.roi;
            if size(val,2) == 4, val = ieRoi2Locs(val); end
        end
    case {'roirect'}
        % sensorGet(sensor,'roi rect')
        % Return ROI as a rect
        if checkfields(sensor,'roi')
            % The data can be stored as a rect or as roiLocs.
            val = sensor.roi;
            if size(val,2) ~= 4, val =  vcLocs2Rect(val); end
        end
    case {'roivolts','roidata','roidatav','roidatavolts'}
        if checkfields(sensor,'roi')
            roiLocs = sensorGet(sensor,'roi locs');
            val = vcGetROIData(sensor,roiLocs,'volts');
        else
            warning('ISET:nosensorroi', 'No roi field. empty returned.');
        end
    case {'roielectrons','roidatae','roidataelectrons'}
        if checkfields(sensor,'roi')
            roiLocs = sensorGet(sensor,'roi locs');
            val = vcGetROIData(sensor,roiLocs,'electrons');
        else warning('ISET:nosensorroi','No roi field. empty returned');
        end
    case {'roivoltsmean'}
        % sensorGet(sensor,'roi volts mean')
        % Mean value for each of the sensor types
        % sensorGet(sensor,'roi volts mean');
        d = sensorGet(sensor,'roi volts');
        if isempty(d), return;
        else 
            nSensor = sensorGet(sensor,'n sensor');
            val = zeros(nSensor,1);
            for ii=1:nSensor
                thisD = d(:,ii);
                val(ii) = mean(thisD(~isnan(thisD)));
            end
        end
    case {'roielectronsmean'}
        % sensorGet(sensor,'roi electrons mean')
        %   Mean value for each of the sensor types
        % sensorGet(sensor,'roi electrons mean');
        % Mean value for each of the sensor types
        % sensorGet(sensor,'roi volts mean');
        d = sensorGet(sensor,'roi electrons');
        if isempty(d), return;
        else 
            nSensor = sensorGet(sensor,'n sensor');
            val = zeros(nSensor,1);
            for ii=1:nSensor
                thisD = d(:,ii);
                val(ii) = mean(thisD(~isnan(thisD)));
            end
        end
    case {'hlinevolts','hlineelectrons','vlinevolts','vlineelectrons'}
        % sensorGet(sensor,'hline volts',row)
        % Returns: val.data and val.pos
        % Each sensor with values on this row in data
        % The positions of the data in pos.
        if isempty(varargin), error('Specify row or col.'); 
        else rc = varargin{1};  % Could be a row or col
        end
        nSensors = sensorGet(sensor,'n sensors');
        
        % Check if the data are in sensor
        if     strfind(param,'volts')
            d = sensorGet(sensor,'volts');
        elseif strfind(param,'electrons')
            d = sensorGet(sensor,'electrons');
        end
        if isempty(d)
            warning('sensorGet:Nolinedata','No data'); 
            return;
        end
        
        support = sensorGet(sensor,'spatial support');
        d = plane2rgb(d,sensor);
        if isequal(param(1),'h')     
            pos = support.x;
        elseif isequal(param(1),'v')
            % To handle 'h' and 'v' case, we transpose the  'v' data to the
            % 'h' format, and we get the y-positions.
            pos = support.y;
            d = imageTranspose(d);
        else error('Unknown orientation.');
        end
        
        % Go get 'em
        val.data   = cell(nSensors,1);
        val.pixPos = cell(nSensors,1);
        for ii=1:nSensors
            thisD = d(rc,:,ii);   % OK because we transposed
            l = find(~isnan(thisD));
            if ~isempty(l)
                val.data{ii} = thisD(l);
                val.pos{ii}  = pos(l)'; 
            end
        end

        % Quantization structure
    case {'quantization','quantizationstructure'}
        val = sensor.quantization;
    case {'nbits','bits'}
        if checkfields(sensor,'quantization','bits')
            val = sensor.quantization.bits;
        end
    case {'max','maxoutput'}
        nbits = sensorGet(sensor,'nbits');
        if isempty(nbits),
            pixel = sensorGet(sensor,'pixel');
            val = pixelGet(pixel,'voltageswing');
        else val = 2^nbits;
        end
    case {'lut','quantizationlut'}
        if checkfields(sensor,'quantization','lut')
            val = sensor.quantization.lut;
        end
    case {'qMethod','quantizationmethod'}
        if checkfields(sensor,'quantization','method')
            val = sensor.quantization.method;
        end

        % Color structure
    case 'color'
        val = sensor.color;
    case {'filterspectra','colorfilters'}
        if sensorCheckHuman(sensor)
            val = sensorGet(sensor,'spectral qe');
        else
            val = sensor.color.filterSpectra;
        end
        
    case {'filternames'}
        % 
        val = sensor.color.filterNames;
    case {'filtercolorletters'}
        % The color letters returned here are in the order of the filter
        % column position in the matrix of filterSpectra, or in the case of
        % human the order of the cone pigments. 
        % 
        % Only the first letter of the filter name is returned.  This
        % information is used in combination with sensorColorOrder to
        % determine plot colors. The letters are a string.
        %
        % The pattern field(see below) describes the position for each
        % filter in the block pattern of color filters.
        names = sensorGet(sensor,'filternames');
        val = zeros(length(names), 1);
        for ii=1:length(names), val(ii) = names{ii}(1); end
        val = char(val);
    case {'filtercolorletterscell'}
        cNames = sensorGet(sensor,'filterColorLetters');
        nFilters = length(cNames);
        val = cell(nFilters,1);
        for ii=1:length(cNames), val{ii} = cNames(ii); end

    case {'filternamescellarray', 'filtercolornamescellarray', ...
            'filternamescell'}
        % N.B.  The order of filter colors returned here corresponds to
        % their position in the columns of filterspectra.  The values in
        % pattern (see below) describes their position in array.
        names = sensorGet(sensor,'filternames');
        val = cell(length(names), 1);
        for ii=1:length(names), val{ii} = char(names{ii}(1)); end
    case {'filterplotcolor','filterplotcolors'}
        % Return an allowable plotting color for this filter, based on the
        % first letter of the filter name.
        % letter = sensorGet(sensor,'filterPlotColor');
        letters = sensorGet(sensor,'filterColorLetters');
        if isempty(varargin), val = letters;
        else                  val = letters(varargin{1});
        end
        % Only return an allowable color.  We could allow w (white) but we
        % don't for now.
        for ii=1:length(val)
            if ~ismember(val(ii),'rgbcmyk'), val(ii) = 'k'; end
        end
    case {'ncolors','nfilters','nsensors','nsensor'}
        val = size(sensorGet(sensor,'filterSpectra'),2);
    case {'ir','infraredfilter','irfilter','otherfilter'}
        % We sometimes put other filters, such as macular pigment, in this
        % slot.  Perhaps we should have an other filter slot.
        if checkfields(sensor,'color','irFilter')
            val = sensor.color.irFilter;
        end
    case {'spectralqe','sensorqe','sensorspectralqe'}
        val = sensorSpectralQE(sensor);

        % There should only be a spectrum associated with the sensor, not
        % with the pixel.  I am not sure how to change over to a single
        % spectral representation, though.  If pixels never existed without
        % an sensor, ... well I am not sure how to get the sensor if only
        % the pixel is passed in.  I am not sure how to enforce
        % consistency. -- BW
    case {'spectrum','sensorspectrum'}
        val = sensor.spectrum;
    case {'wave','wavelength'}
        val = sensor.spectrum.wave(:);
    case {'binwidth','waveresolution','wavelengthresolution'}
        wave = sensorGet(sensor,'wave');
        if length(wave) > 1, val = wave(2) - wave(1);
        else val = 1;
        end
    case {'nwave','nwaves','numberofwavelengthsamples'}
        val = length(sensorGet(sensor,'wave'));

        % Color filter array quantities
    case {'cfa','colorfilterarray'}
        val = sensor.cfa;

        % I removed the unitBlock data structure because everything that
        % was in unit block can be derived from the cfa.pattern entry.  We
        % are coding the cfa.pattern entry as a small matrix.  So, for
        % example, if it is a 2x2 Bayer pattern, cfa.pattern = [1 2; 2 3]
        % for a red, green, green, blue pattern.  The former entries in
        % unitBlock are redundant with this value and the pixel size.  So,
        % we got rid of them.
    case {'unitblockrows'}
        % sensorGet(sensor,'unit block rows')
        
        % Human patterns don't have block sizes.
        if sensorCheckHuman(sensor), val=1;
        else val = size(sensorGet(sensor,'pattern'),1);
        end
        
    case 'unitblockcols'
        % sensorGet(sensor,'unit block cols')
        
        % Human patterns don't have block sizes.
        if sensorCheckHuman(sensor), val=1;
        else val = size(sensorGet(sensor,'pattern'),2);
        end
        
    case {'cfasize','unitblocksize'}
        % We use this to make sure the sensor size is an even multiple of
        % the cfa size. This could be a pair of calls to cols and rows
        % (above).
        
        % Human patterns don't have block sizes.
        if sensorCheckHuman(sensor), val= [1 1];
        else    val = size(sensorGet(sensor,'pattern'));
        end
        
    case 'unitblockconfig'
        % val = sensor.cfa.unitBlock.config;
        % Is this still used?
        pixel = sensorGet(sensor,'pixel');
        p = pixelGet(pixel,'pixelSize','m');
        [X,Y] = meshgrid((0:(size(cfa.pattern,2)-1))*p(2), ...
                         (0:(size(cfa.pattern,1)-1))*p(1));
        val = [X(:),Y(:)];
    
    case {'patterncolors','pcolors','blockcolors'}
        % patternColors = sensorGet(sensor,'patternColors');
        % Returns letters suggesting the color of each pixel
        
        pattern = sensorGet(sensor,'pattern');  %CFA block
        filterColorLetters = sensorGet(sensor,'filterColorLetters');
        knownColorLetters = sensorColorOrder('string');
        knownFilters = ismember(filterColorLetters,knownColorLetters);
        % Assign unknown color filter strings to black (k).
        l = find(~knownFilters, 1);
        if ~isempty(l), filterColorLetters(l) = 'k'; end
        % Create a block that has letters instead of numbers
        val = filterColorLetters(pattern);

    case {'pattern','cfapattern'}
        if checkfields(sensor,'cfa','pattern')
            val = sensor.cfa.pattern;
        end
    case 'cfaname'
        % We look up various standard names
        val = sensorCFAName(sensor);

        % Pixel related parameters
    case 'pixel'
        val = sensor.pixel;
        
    case {'dr','dynamicrange','sensordynamicrange'}
        val = sensorDR(sensor);

    case 'diffusionmtf'
        val = sensor.diffusionMTF;

        % Pixel-wise FPN parameters
    case {'fpnparameters','fpn','fpnoffsetgain','fpnoffsetandgain'}
        val = [sensorGet(sensor,'sigmaOffsetFPN'), ...
               sensorGet(sensor,'sigmaGainFPN')];
    case {'dsnulevel','sigmaoffsetfpn','offsetfpn', ...
          'offset','offsetsd','dsnusigma','sigmadsnu'}
        % This value is stored in volts
        val = sensor.sigmaOffsetFPN;
    case {'sigmagainfpn', 'gainfpn', 'gain', 'gainsd', ...
          'prnusigma', 'sigmaprnu', 'prnulevel'}
        % This is a percentage, between 0 and 100, always.
        val = sensor.sigmaGainFPN;
    
    case {'dsnuimage','offsetfpnimage'} % Dark signal non uniformity (DSNU)
        % These should probably go away because we compute them afresh
        % every time.
        if checkfields(sensor,'offsetFPNimage')
            val = sensor.offsetFPNimage;
        end
    case {'prnuimage','gainfpnimage'}  % Photo response non uniformity
        % These should probably go away because we compute them afresh
        % every time.
        if checkfields(sensor,'gainFPNimage')
            val = sensor.gainFPNimage;
        end

        % These are column-wise FPN parameters
    case {'columnfpn','columnfixedpatternnoise','colfpn'}
        % This is stored as a vector (offset,gain) standard deviations in
        % volts.  This is unlike the storage format for array dsnu and prnu
        if checkfields(sensor,'columnFPN'), val = sensor.columnFPN; 
        else
            val = [0,0];
        end
    case {'columndsnu','columnfpnoffset','colfpnoffset','coldsnu'}
        tmp = sensorGet(sensor,'columnfpn'); val = tmp(1);
    case {'columnprnu','columnfpngain','colfpngain','colprnu'}
        tmp = sensorGet(sensor,'columnfpn'); val = tmp(2);
    case {'coloffsetfpnvector','coloffsetfpn','coloffset'}
        if checkfields(sensor,'colOffset'), val = sensor.colOffset; end
    case {'colgainfpnvector','colgainfpn','colgain'}
        if checkfields(sensor,'colGain'),val = sensor.colGain; end
        
        % Noise management
    case {'noiseflag','shotnoiseflag'}
        % 0 means no noise
        % 1 means shot noise but no electronics noise
        % 2 means shot noise and electronics noise
        if checkfields(sensor,'noiseFlag'), val = sensor.noiseFlag;
        else val = 2;    % Compute both electronic and shot noise
        end
    case {'reusenoise'}
        if checkfields(sensor,'reuseNoise'), val = sensor.reuseNoise;
        else val = 0;    % Do not reuse
        end
    case {'noiseseed'}
        if checkfields(sensor,'noiseSeed'), val = sensor.noiseSeed; end

    case {'ngridsamples', 'pixelsamples', ...
          'nsamplesperpixel', 'npixelsamplesforcomputing'}
        % Default is 1. If parameter is not set, we return the default.
        if checkfields(sensor,'samplesPerPixel')
            val = sensor.samplesPerPixel;
        else val = 1;
        end

        % Exposure related
    case {'exposuremethod','expmethod'}
        % We plan to re-write the exposure parameters into a sub-structure
        % that lives inside the sensor, sensor.exposure.XXX
        tmp = sensorGet(sensor,'exptimes');
        p   = sensorGet(sensor,'pattern');
        if     isscalar(tmp), val = 'singleExposure';
        elseif isvector(tmp),  val = 'bracketedExposure';
        elseif isequal(size(p),size(tmp)),  val = 'cfaExposure';
        end
    case {'integrationtime','integrationtimes','exptime','exptimes', ...
          'exposuretimes','exposuretime','exposureduration', ...
          'exposuredurations'}
        % This can be a single number, a vector, or a matrix that matches
        % the size of the pattern slot. Each one of these cases is handled
        % differently by sensorComputeImage.  The units are seconds by
        % default.
        % sensorGet(sensor,'expTime','s')
        % sensorGet(sensor,'expTime','us')
        val = sensor.integrationTime;
        if ~isempty(varargin)
            val = val*ieUnitScaleFactor(varargin{1});
        end
    case {'uniqueintegrationtimes','uniqueexptime','uniqueexptimes'}
        val = unique(sensor.integrationTime);
    case {'centralexposure','geometricmeanexposuretime'}
        % We return the geometric mean of the exposure times
        % We should consider making this the geometric mean of the unique
        % exposures.
        eTimes = sensorGet(sensor,'exptimes');
        val = prod(eTimes(:))^(1/length(eTimes(:)));
    case {'autoexp','autoexposure','automaticexposure'}
        val = sensor.AE;
    case {'nexposures'}
        % We can handle multiple exposure times.
        val = numel(sensorGet(sensor,'expTime'));
    case {'exposureplane'}
        % When there are multiple exposures, show the middle integration
        % time, much like a bracketing idea.
        % N.B. In the case where there is a different exposure for every
        % position in the CFA, we wouldn't normally use this.  In that case
        % we only have a single integrated CFA.
        if checkfields(sensor,'exposurePlane'), val = sensor.exposurePlane;
        else val = floor(sensorGet(sensor,'nExposures')/2) + 1;
        end

    case {'cds','correlateddoublesampling'}
        val = sensor.CDS;

        % Microlens related
    case {'vignettingflag','vignetting', ...
          'pixelvignetting','bareetendue', ...
          'sensorbareetendue','nomicrolensetendue'}
        % If the vignetting flag has not been set, treat it as 'skip',
        % which is 0.
        if checkfields(sensor,'data','vignetting'),
            if isempty(sensor.data.vignetting), val = 0;
            else                             val = sensor.data.vignetting;
            end
        else
            val = 0;
        end
    case {'vignettingname'}
        pvFlag = sensorGet(sensor,'vignettingFlag');
        switch pvFlag
            case 0
                val = 'skip';
            case 1
                val = 'bare';
            case 2
                val = 'centered';
            case 3
                val = 'optimal';
            otherwise
                error('Bad pixel vignetting flag')
        end

    case {'microlens','ulens','mlens','ml'}
        if checkfields(sensor,'ml'), val = sensor.ml; end

        % Field of view and sampling density
    case {'fov','sensorfov','fovhorizontal','fovh','hfov'}
        % sensorGet(sensor,'fov',sDist,oi); - Explicit scene dist in m
        % sensorGet(sensor,'fov',scene,oi); - Explicit scene
        % sensorGet(sensor,'fov');          - Uses defaults.  Dangerous.
        %
        % This is the horizontal field of view (default)
        % We compute it from the distance between the lens and the sensor
        % surface and we also use the sensor array width.
        % The assumption here is that the sensor is at the proper focal
        % distance for the scene.  If the scene is at infinity, then the
        % focal distance is the focal length.  But if the scene is close,
        % then we might correct.
        %
        if ~isempty(varargin), scene = varargin{1};
        else                   scene = vcGetObject('scene');
        end
        if length(varargin) > 1, oi = varargin{2};
        else                     oi = vcGetObject('oi');
        end
        
        % If no scene is sent in, assume the scene is infinitely far away.
        if isempty(scene), sDist = Inf;
        elseif isstruct(scene), sDist = sceneGet(scene,'distance');
        else   sDist = scene; % sDist is directly passed in (in meters)
        end
        
        % If there is no oi, then use the default optics focal length. The
        % image distance depends on the scene distance and focal length via
        % the lensmaker's formula, (we assume the sensor is at the proper
        % focal distance).
        if ~isempty(oi) 
            distance = oiGet(oi,'optics image distance', sDist);
        elseif sensorCheckHuman(sensor)
            distance = oiGet(oiCreate('human'), 'focal length');
        else
            distance = opticsGet(opticsCreate,'focal length');
            fprintf('Fov estimated using focal length = %f m\n',distance);
        end
        width = sensorGet(sensor, 'arraywidth');
        val = 2 * atand(0.5 * width / distance);
    case {'fovvertical','vfov','fovv'}
        % This is  the vertical field of view
        if ~isempty(varargin), scene = varargin{1};
        else                   scene = vcGetObject('scene');
        end
        if length(varargin) > 1, oi = varargin{2};
        else                     oi = vcGetObject('oi');
        end
        
        % If no scene is sent in, assume the scene is infinitely far away.
        if isempty(scene), sDist = Inf;
        elseif isstruct(scene), sDist = sceneGet(scene,'distance');
        else   sDist = scene; % sDist is directly passed in (in meters)
        end
        
        % If there is no oi, then use the default optics focal length. The
        % image distance depends on the scene distance and focal length via
        % the lensmaker's formula, (we assume the sensor is at the proper
        % focal distance).
        if ~isempty(oi) 
            distance = oiGet(oi,'optics image distance', sDist);
        elseif sensorCheckHuman(sensor)
            distance = oiGet(oiCreate('human'), 'focal length');
        else
            distance = opticsGet(opticsCreate,'focal length');
            fprintf('Fov estimated using focal length = %f m\n',distance);
        end
        height = sensorGet(sensor, 'arrayheight');
        val = 2*atand(0.5 * height / distance);
    case {'hdegperpixel','degpersample','degreesperpixel'}
        % degPerPixel = sensorGet(sensor,'h deg per pixel',oi);
        % Horizontal field of view divided by number of pixels
        sz =  sensorGet(sensor,'size');
        
        if isempty(varargin), oi = vcGetObject('oi');
        else oi = varargin{1};
        end
        
        % The horizontal field of view should incorporate information from
        % the optics.
        val = sensorGet(sensor,'hfov',oi)/sz(2);
    case {'vdegperpixel','vdegreesperpixel'}
        sz =  sensorGet(sensor,'size');
        val = sensorGet(sensor,'vfov')/sz(1);
    case {'hdegperdistance','degperdistance'}
        % sensorGet(sensor,'h deg per distance','mm')
        % sensorGet(sensor,'h deg per distance','mm',scene,oi);
        % Degrees of visual angle per meter or other spatial unit
        if isempty(varargin), unit = 'm'; else unit = varargin{1}; end
        width = sensorGet(sensor,'width',unit);
        
        if length(varargin) < 2, scene = vcGetObject('scene');
        else scene = varargin{2};
        end
        
        % We want the optics to do this right.
        if length(varargin) < 3, oi = vcGetObject('oi');
        else oi = varargin{3};
        end
        
        fov   =  sensorGet(sensor,'fov',scene, oi);
        val   = fov/width;
        
    case {'vdegperdistance'}
        % sensorGet(sensor,'v deg per distance','mm') Degrees of visual
        % angle per meter or other spatial unit
        if isempty(varargin), unit = 'm'; else unit = varargin{1}; end
        width = sensorGet(sensor,'height',unit);
        fov =  sensorGet(sensor,'vfov');
        val = fov/width;
        
        % Computational flags
    case {'sensorcompute','sensorcomputemethod'}
        % Swap in a sensorCompute routine.  If this is empty, then the
        % standard vcamera\sensor\mySensorCompute routine will be used.
        if checkfields(sensor,'sensorComputeMethod')
            val = sensor.sensorComputeMethod;
        else  val = 'mySensorCompute';
        end
    case {'consistency','computationalconsistency'}
        % If the consistency field is not present, assume false and set it
        % false.  This checks whether the parameters and the displayed
        % image are consistent/updated.
        if checkfields(sensor,'consistency'), val = sensor.consistency;
        else sensorSet(sensor,'consistency',0); val = 0;
        end

    case {'mccrecthandles'}
        % These are handles to the squares on the MCC selection regions
        % see macbethSelect
        if checkfields(sensor,'mccRectHandles')
            val = sensor.mccRectHandles;
        end
    case {'mccpointlocs','mcccornerpoints'}
        % Corner points for the whole MCC chart
        if checkfields(sensor,'mccCornerPoints')
            val = sensor.mccCornerPoints;
        end

        % Display image
    case {'rgb', 'rgbimage'}
        % sensorGet(sensor,'rgb',dataType,gam,scaleMax)
        dataType = 'volts';
        gam = 1;
        scaleMax = 0;
        if ~isempty(varargin), dataType = varargin{1}; end
        if length(varargin) > 1, gam = varargin{2}; end
        if length(varargin) > 2, scaleMax = varargin{3}; end
        val = sensorData2Image(sensor,dataType,gam,scaleMax);

        % Human sensor (eye) case
    case {'human'}
        % Structure containing information about human cone case
        % Only applies when the name field has the string 'human' in it.
        if checkfields(sensor,'human'), val = sensor.human; end
    case {'humancone'}
        if checkfields(sensor,'human')
            val = sensor.human.cone;
        end
    case {'humanlens', 'lens'}
        if checkfields(sensor,'human', 'lens')
            val = sensor.human.lens;
        end
    case {'humanlenstrans', 'humanlenstransmittance', 'lenstransmittance'}
        % Lens transmittance
        lens = sensorGet(sensor, 'human lens');
        val = lensGet(lens, ' transmittance');
    case {'humanlensabsorption','lensabsorption','lensabsorptance'}
        lens = sensorGet(sensor, 'human lens');
        val = lensGet(lens, 'absorption');
        
    case {'humanmaculardensity', 'macdens','maculardensity'}
        macular = sensorGet(sensor, 'human macular');
        val = macularGet(macular, 'density');
    case {'humanmaculartransmittance', 'maculartrans', ...
          'maculartransmittance'}
        % senosrGet(sensor,'macular transmittance')
        macular = sensorGet(sensor, 'human macular');
        val = macularGet(macular, 'transmittance');
    case {'humanmacularabsorption', 'macularabsorption'}
        % sensorGet(sensor,'macular absorption')
        macular = sensorGet(sensor, 'human macular');
        val = macularGet(macular, 'absorption');
    case {'humanmacular', 'macular'}
        if checkfields(sensor,'human', 'macular')
            val = sensor.human.macular;
        end
    case {'humanoculartransmittance'}
        lens = sensorGet(sensor, 'human lens');
        macular = sensorGet(sensor, 'human macular');
        val = lensGet(lens,'transmittance') .* ...
                macularGet(macular, 'transmittance');
    case {'humaneffectiveabsorptance'}
        % Combines cone photopigment, ocular (lens and macular)
        % transmittance and peak efficiency
        cone        = sensorGet(sensor, 'human cone');
        absorptance = coneGet(cone, 'absorptance');
        eyeTrans    = sensorGet(sensor, 'human ocular transmittance');
        peakEfficiency = coneGet(cone,'peak efficiency');
        
        val = (absorptance.*repmat(eyeTrans, [1 size(absorptance,2)]))* ...
            diag(peakEfficiency);
        nFilters = length(sensorGet(sensor, 'filter names'));
        if size(val, 2) ~= nFilters
            val = [zeros(size(val,1),1) val];
        end
        assert(size(val, 2) == nFilters, 'size mismatch: filter name/qe');
    case {'humanconetype','conetype'}
        % Blank (K) K=1 and L,M,S cone at each position
        % L=2, M=3 or S=4 (K means none)
        % Some number of cone types as cone positions.
        if checkfields(sensor,'human','coneType')
            val = sensor.human.coneType;
        end
    case {'humanconedensities','conedensity', 'densities'}
        %- densities used to generate mosaic (K,L,M,S)
        val = coneGet(sensorGet(sensor, 'human cone'), 'spatial density');
    case   {'humanconelocs','conexy','conelocs','xy'}
        %- xy position of the cones in the mosaic
        if checkfields(sensor,'human','xy'), val = sensor.human.xy; end
    case {'humanrseed','humanconeseed', 'rseed'}
        % random seed for generating cone mosaic
        % Should get rid of humanrseed alias
        if checkfields(sensor,'human','rSeed')
            val = sensor.human.rSeed;
        end
        
    case {'sampletimeinterval', 'timeinterval'}
        if checkfields(sensor, 'human', 'timeInterval')
            val = sensor.human.timeInterval;
        end
        
    case {'adaptationgain'}
        % Adaptation gain
        if checkfields(sensor,'human','adaptGain')
            val = sensor.human.adaptGain;
        end
    case {'adaptationoffset'}
        % Adaptation offset
        if checkfields(sensor, 'human', 'adaptOffset')
            val = sensor.humna.adaptOffset;
        end
        
    case {'adapteddata'}
        % volts after cone adaptation
        gainMap = sensorGet(sensor, 'adaptation gain');
        offset  = sensorGet(sensor, 'adaptation volts');
        volts   = sensorGet(sensor, 'volts');
        
        if isscaler(gainMap)
            val = volts .* gainMap - offset;
        elseif ismatrix(gainMap)
            nSamples = size(volts, 3); 
            val = volts .* repmat(gainMap, [1 1 nSamples]) - offset;
        else
            val = volts .* gainMap - offset;
            val(isnan(val)) = 0;
        end

    case {'movementpositions', 'sensorpositions', 'positions'}
        % Nx2 vector of (x,y) positions in deg
        if checkfields(sensor,'movement','pos'), val = sensor.movement.pos;
        else val = [0,0]; 
        end
    case {'positionsx','sensorpositionsx'}
        if checkfields(sensor,'movement','pos')
            val = sensor.movement.pos(:,1);
        else val = 0;
        end
    case {'positionsy','sensorpositionsy'}
        if checkfields(sensor,'movement','pos')
            val = sensor.movement.pos(:,2);
        else val = 0;
        end
    case {'framesperposition','exposuretimesperposition','etimeperpos'}
        % Exposure frames for each (x,y) position
        % This is a vector with some number of exposures for each x,y
        % position (deg)
        warning('This field might be removed in the future');
        disp(['For a general eyemovement sequence '...
            'set sensor positions to the sensor.']);
        
        if checkfields(sensor,'movement','framesPerPosition')
            val = sensor.movement.framesPerPosition;
        else
            val = 1;
        end
        
    case {'tottime', 'totaltime'}
        % total time of exposure
        sampTime = sensorGet(sensor, 'time interval');
        seqLen = size(sensorGet(sensor, 'sensor positions'), 1);
        val = seqLen * sampTime;
        
    case {'eyemove', 'eyemovement'}
        % eye movement structure
        if checkfields(sensor, 'human', 'eyemove')
            val = sensor.human.eyemove;
        end
    case {'emflag'}
        % eye movement type
        if checkfields(sensor, 'human', 'eyemove', 'emFlag')
            val = sensor.human.eyemove.emFlag;
        end
    case {'emtremor'}
        % eye movemenet tremor structure
        if checkfields(sensor, 'human', 'eyemove', 'tremor')
            val = sensor.human.eyemove.tremor;
        end        
    case {'emdrift'}
        % eye movement drift structure
        if checkfields(sensor, 'human', 'eyemove', 'drift')
            val = sensor.human.eyemove.drift;
        end
    case {'emmsaccade', 'emmicrosaccade'}
        % eye movement microsaccade structure
        if checkfields(sensor, 'human', 'eyemove', 'msaccade')
            val = sensor.human.eyemove.msaccade;
        end
    otherwise
        error('Unknown sensor parameter.');
end

end

%--------------------------------
function cfaName = sensorCFAName(sensor)
% Determine the cfa name in order to populate the lists in pop up boxes.
%
%     cfaName = sensorCFAName(sensor)
%
% If sensor is passed in, return a standard name for the CFA types.
% If sensor is empty or absent return the list of standard names.
%
% The normal naming convention for CFA is to read left to right.  For
% example,
%     G B
%     R G
% is coded as 'gbrg'
% The pattern matrix stored in sensor is a 2x2 array, usually.  Thus, the
% values are stored as
%     [2 1; 3 2]
%   =    2 1
%        3 2
%
% This is a bit confusing.  Sorry.
%
% Examples:
%   cfaName = sensorCFAName(sensor)
%   cfaNames = sensorCFAName
%

if notDefined('sensor')
    cfaName = sensorCFANameList;
    return;
end

p = sensorGet(sensor,'pattern');
filterColors = sensorGet(sensor,'filterColorLetters');
filterColors = sort(filterColors);

if length(p) == 1
    cfaName = 'Monochrome';
elseif length(p) > 4
    cfaName = 'Other'; return;
elseif isequal(p,[ 2 1; 3 2]) && strcmp(filterColors,'bgr');
    cfaName = 'bayer-grbg';
elseif isequal(p,[ 1 2; 2 3]) && strcmp(filterColors,'bgr');
    cfaName = 'bayer-rggb';
elseif isequal(p,[ 2 3; 1 2]) && strcmp(filterColors,'bgr');
    cfaName = 'bayer-gbrg';
elseif isequal(p,[ 3 2; 2 1]) && strcmp(filterColors,'bgr');
    cfaName = 'bayer-bggr';
elseif isequal(p,[ 2 3; 1 2]) && strcmp(filterColors,'cym');
    cfaName = 'bayer-ycmy';
elseif  isequal(p,[ 1 3; 2 4])
    cfaName = 'Four Color';
else
    cfaName = 'Other';
end

end

%------------------------------------------
function spectralQE = sensorSpectralQE(sensor)
% Compute the sensor spectral QE
%
%    spectralQE = sensorSpectralQE(sensor)
%
% Combine the pixel detector, the sensor color filters, and the infrared
% filter into a sensor spectral QE.   If the variable wave is in the
% calling arguments, the spectralQE is returned interpolated to the
% wavelengths in wave.
%

if sensorCheckHuman(sensor)
    % This combines the lens, macular pigment and cone pigments into a
    % single, overall spectral QE.
    spectralQE = sensorGet(sensor, 'human effective absorptance');

else
    
    sensorIR = sensorGet(sensor,'ir filter');
    cf       = sensorGet(sensor,'filter spectra');
    
    pixelQE = pixelGet(sensor.pixel,'qe');
    if isempty(pixelQE)
        warndlg('Empty pixel QE. Assuming QE(lambda) = 1.0');
        pixelQE = ones(size(sensorIR(:)));
    end
    
    % Compute the combined wavelength sensitivity including the ir filter,
    % the pixel QE, and the color filters.
    spectralQE = diag(pixelQE(:) .* sensorIR(:)) * cf;
end

end

%------------------------
function val = sensorColorData(data,sensor,whichSensor)
% Retrieve data from one of the sensor planes.
%
% The data are returned in a vector, not a plane.
%
% This should also be able to return a planar image with zeros, not just a
% vector.  But surely this form was written some time ago.

% In most cases, we can convert the data to a 3D and return the values in
% the RGB 3D.  In human, we have the whole pattern.  Probably we should
% always get the
%
% This might work in both cases ... but sensorDetermineCFA may not be
% running right for the human case.  Have a look.  For now, we are only
% using the 'ideal' condition with human.
%
% electrons = sensorGet(sensor,'electrons');
% [tmp,n] = sensorDetermineCFA(sensor);
% b = (n==whichSensor);
% val = electrons(b);

% The one we have been using
rgb        = plane2rgb(data,sensor);
thisSensor = rgb(:,:,whichSensor);
l   = ~isnan(thisSensor);
val = thisSensor(l);

end

% TODO:
%
% The pixel height and width may differ from the center-to-center distance
% between pixels in the sensor.  This is because of additional gaps between
% the pixels (say for wires). The center-to-center spacing between pixels
% is contained in the deltaX, deltaY parameters.
%
% For consistency with the other structures, I also introduced the
% parameters hspatialresolution and wspatialresolution to refer to the
% deltaY and deltaX center-to-center spacing of the pixels.  We might
% consider adding just hresolution and wresolution for spatial, with
% angular being special.
%
% In several cases we use the spatial coordinates of the sensor array to
% define a coordinate system, say for the optical image.  In this case, the
% integer coordinates are defined by the deltaX and deltaY values.
%
% get cfa matrix as letters or numbers via sensorDetermineCFA in here.

