function val = pixelGet(pixel,param,varargin)
%Get data from the pixel structure
%
%     val = pixelGet(pixel,param,varargin)
%
% The list of properties is below.  For many of the size properties, you
% can request different spatial units by specifying the desired unit in
% a calling argument.  For example, pixelGet(pixel,'pixelWidth','um')
% returns the pixel width in microns.  
%
% For this and other get routines, capitalization is irrelevant.
%
% The get routines usually have several aliases for the same quantity.
% Read the code to see the aliases.
%
% Pixel spatial size
%      {'pixelwidth'}          - pixel width (meters)
%      {'pixelheight'}         - pixel height (meters)
%      {'pixelwidthgap'}       - width gap between pixels (meters) 
%      {'pixelheightgap'}      - height gap between pixels (meters)
%      {'pixelsize'}           - (width,height) vector
%      {'pixelarea'}           - width*height (meters^2)
%      {'wspatialresolution'}  - spacing in x-direction (width spatial resolution)    
%      {'hspatialresolution'}  - spacing in y-direction (height spatial resolution)
%      {'xyspacing'}           - dimension is (x,y) but size is (row,col)
%
% Photodetector properties
%      {'pdsize'}              -  size (height,width)
%      {'photodetectorwidth'}  -  width (meters)
%      {'photodetectorheight'} -  height (meters)
%      {'fillfactor'}          -  fraction of pixel area occupied by photodetector
%      {'pixeldepth'}          -  depth (meters)
%      {'pdxpos'}              -  x-position inside pixel (meters)
%      {'pdypos'}              -  x-position inside pixel (meters)
%      {'pdposition'}          -  (x,y)-position inside pixel (meters)
%      {'pddimension'}         -  size in (x,y) format (size is (row,col) format)
%      {'pdarea'}              -  area (m^2)
%
% Optical properties
%      {'refractiveindices'}   - refractive indices of air, materials, silicon
%      {'layerthicknesses'}    - thickness of different materials (meters)
%      {'stackheight'}         - thickness of different materials (meters)
%      {'pixelspectrum'}       - spectrum structure for pixel
%      {'wavelength'}          - wavelength sampling for pixel 
%      {'binwidth'}            - bin width of wavelength 
%      {'nwave'}               - number of wavelength samples
%      {'pdspectralqe'}        - photodetector spectral quantum efficiency
%
% Electrical properties
%      {'conversiongain'}      - volts per electron 
%      {'voltageswing'}        - maximum voltage 
%      {'wellcapacity'}        - maximum number of electrons (=vSwing/convGain)
%      {'darkcurrentdensity'}  - dark current in (amps/m^2)
%      {'darkcurrent'}         - dark current per photodetector (=dkDensity*pdarea)
%      {'darkvoltage'}         - dark voltage per photodetector
%      {'darkelectrons'}       - electrons / pixel / second
%      {'readnoiseelectrons'}  - standard deviation of read noise in electrons
%      {'readnoisevolts'}      - standard deviation of read noise volts
%      {'readnoisemillivolts'} - standard deviation of read noise in millivots
%      {'pdspectralsr'}        - photodetector spectral responsivity(pixelSR)
%      {'pixeldr'}             - pixel dynamic range.
%
% Examples:
%   width = pixelGet(pixel,'width','um')
%   sz = pixelGet(pixel,'size');
%   dkV = pixelGet(pixel,'darkVoltage')
%  
% Copyright ImagEval Consultants, LLC, 2005.

if ~exist('pixel','var') || isempty(pixel), error('Must define pixel.'); end
if ~exist('param','var') || isempty(param), error('Must define parameter.'); end

param = ieParamFormat(param);

switch param
    case {'width','pixelwidth','pixelwidthmeters'}    %M
        % pixelGet(pixel,'width','microns')
        val = pixel.width;
        if ~isempty(varargin), val = val*ieUnitScaleFactor(varargin{1}); end
    case {'height','pixelheight','pixelheightmeters'}   %M
        val = pixel.height;
        if ~isempty(varargin), val = val*ieUnitScaleFactor(varargin{1}); end
    case {'pixeldepth','depth','pixeldepthmeters'}   %M
        val = sum(pixelGet(pixel,'layerthickness'));
        if ~isempty(varargin), val = val*ieUnitScaleFactor(varargin{1}); end
        
    case {'pixelwidthgap','widthgap'} %M
        val = pixel.widthGap;
        if ~isempty(varargin), val = val*ieUnitScaleFactor(varargin{1}); end
        
    case {'pixelheightgap','heightgap'}    %M
        val = pixel.heightGap;
        if ~isempty(varargin), val = val*ieUnitScaleFactor(varargin{1}); end
        
    case {'wspatialresolution','deltax'}    
        val = pixelGet(pixel,'width') + pixelGet(pixel,'widthGap');
        if ~isempty(varargin), val = val*ieUnitScaleFactor(varargin{1}); end
        
    case {'hspatialresolution','deltay'}
        val = pixelGet(pixel,'height') + pixelGet(pixel,'heightGap');
        if ~isempty(varargin), val = val*ieUnitScaleFactor(varargin{1}); end
        
    case {'xyspacing','dimension'}    % Note:  dimension is (x,y) but size is (row,col)
        val = [pixelGet(pixel,'deltax'), pixelGet(pixel,'deltay')];
        if ~isempty(varargin), val = val*ieUnitScaleFactor(varargin{1}); end

    case {'pixelsize','size'}
        val = [pixelGet(pixel,'deltay'), pixelGet(pixel,'deltax')];
        if ~isempty(varargin), val = val*ieUnitScaleFactor(varargin{1}); end
        
    case {'pixelarea','area'}
        % This includes the gap dimension.  
        val = prod(pixelGet(pixel,'size'));
        
        % Photodetector sizes and positions
    case {'photodetectorwidth','pdwidth'}  %M
        val = pixel.pdWidth;
        if ~isempty(varargin), val = val*ieUnitScaleFactor(varargin{1}); end
        
    case {'photodetectorheight','pdheight'} %M
        val = pixel.pdHeight;
        if ~isempty(varargin), val = val*ieUnitScaleFactor(varargin{1}); end
        
    case {'pdsize','photodetectorsize'}
        % pixelGet(pixel,'pd height','um')
        % pixelGet(pixel,'pd height','m')
        val = [pixelGet(pixel,'pdHeight'),pixelGet(pixel,'pdWidth')];
        if ~isempty(varargin), val = val*ieUnitScaleFactor(varargin{1}); end
        
    case {'layerthickness','layerthicknesses'} %M
        val = pixel.layerThickness;
        if ~isempty(varargin), val = val*ieUnitScaleFactor(varargin{1}); end
        
    case {'stackheight'} %M
        val = sum(pixel.layerThickness);
        if ~isempty(varargin), val = val*ieUnitScaleFactor(varargin{1}); end 
        
    case {'refractiveindex','refractiveindices','n'} %
        val = pixel.n;
        
        % The photodetector x and y positions are in some units ... please
        % comment here.  Is it in Meters?  I think so.
    case 'pdxpos'       %M?
        val = pixel.pdXpos;
    case 'pdypos'
        val = pixel.pdYpos;
    case 'pdposition'
        val = [pixelGet(pixel,'pdXpos'),pixelGet(pixel,'pdYpos')];
        
    case 'pddimension'  % Dimension is (x,y) and size is (row,col)
        val = [pixelGet(pixel,'pdWidth'),pixelGet(pixel,'pdHeight')];
        if ~isempty(varargin), val = val*ieUnitScaleFactor(varargin{1}); end
        
  
    case 'pdarea'
        val = prod(pixelGet(pixel,'pdsize')); 
    case 'fillfactor';
        val = pixelGet(pixel,'pdArea') / pixelGet(pixel,'area');

        % Electrical properties
    case 'conversiongain'  % Volts/e-
        val = pixel.conversionGain;
    case {'voltageswing','vswing'}     % Volts
        val = pixel.voltageSwing;
    case {'wellcapacity'}   % In electrons
        val = pixelGet(pixel,'voltageSwing')/pixelGet(pixel,'conversionGain');
        
    case {'darkcurrentdensity'}     %Amps/m2
        % Amps/pdetector * (pdetector/M2) = Amps/m2
        val = pixelGet(pixel,'darkcurrent')/pixelGet(pixel,'pdarea');
        
    case {'darkcurrent','darkcurrentperpixel'}            % Amps/photodetector
        % V/sec*(1/(V/e-))*charge/e- = charge/sec at each pixel
        val = (pixelGet(pixel,'darkvoltage')/pixelGet(pixel,'conversiongain'))*vcConstants('q');
        % val = pixelGet(pixel,'darkcurrentdensity') * pixelGet(pixel,'pdarea');
        
    case {'darkvolt','darkvoltage','darkvolts','darkvoltageperpixelpersec'}
        % V/e * (Charge/sec)/M^2 * M^2 * (1/Charge/electron) = V/sec
        %  val = ...
        %   pixelGet(pixel,'conversiongain')*pixelGet(pixel,'darkcurrent')/vcConstants('q');
        val = pixel.darkVoltage;        %Volts per second
        
    case {'darkelectrons'}
        % Electrons / pixel / second
        val = pixelGet(pixel,'darkvoltage')/pixelGet(pixel,'conversiongain');

    case {'readnoiseelectrons'}         %standard deviation in electrons
        val = pixel.readNoise/pixelGet(pixel,'conversiongain');
    case {'readnoisevolts','readnoise'} %standard deviation in Volts
        val = pixel.readNoise;
    case {'readnoisemillivolts'}
        val = pixel.readNoise*10^3;
        
    case {'spectrum','pixelspectrum'}
        val = pixel.spectrum;
    case {'wave','wavelength','wavelengthsamples'}         %nm
        val = pixel.spectrum.wave(:);
    case {'binwidth','wavelengthresolution'}     %nm
        wave = pixelGet(pixel,'wave');
        if length(wave) > 1, val = wave(2) - wave(1);
        else val = 1;
        end
    case {'nwave','nwaves','numberofwavelengthsamples'}
        val = length(pixel.spectrum.wave);
        
    case {'spectralqe','pdspectralqe','qe','pixelspectralqe'}
        val = pixel.spectralQE(:);
    case {'spectralsr','pdspectralsr','sr'}
        % Converts QE to SR
        val = pixelSR(pixel);
        
    case {'pixeldr','pixeldynamicrange'}
        val = pixelDR(ISA);
        
    otherwise
        error('Unknown param: %s',param);
end            
