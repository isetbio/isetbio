function [sensor, coneP] = sensorCreate(sensorName,pixel,varargin)
% Create an image sensor array structure
%
%   [sensor,coneP] = sensorCreate(sensorName,[pixel/coneP],varargin)
%
% The sensor array comprises a matrix of cone types.  The distribution and
% spatial array of cone types is specified by the parameter coneP.  This
% parameter is used in the function sensorCreateConeMosaic.
%
% Several type of image sensors can be created, including multispectral and
% a model of the human cone mosaic.
%
% Human cone mosaic
%      {'human'} - Uses Stockman Quanta LMS cones
%                  Default params:  
%                     params.XXX
%
%  Bayer RGB combinations - may be deprecated
%      {'bayer-grbg'}
%      {'bayer-rggb'}
%      {'bayer-bggr'}
%      {'bayer-gbrg'}
%
% Other types
%      {'monochrome'}
%      {'monochrome array'} - sensorCreate('monochrome array',N);
%
% Multiple channel sensors can be created
%      {'grbc'}        - green, red, blue, cyan
%      {'interleaved'} - One transparent channel and 3 RGB.  Same as RGBC
%                        or RGBW
%
% See also: sensorCreateConeMosaic, coneCreate, sensorConePlot
%
% Examples
%  Basic cone mosaic
%
%   coneP = coneCreate;               % Specify cone properties
%   sensor = sensorCreate('human');   
%   sensor = sensorCreate('human',coneP);
%
%   sensorConePlot(sensor)
%
%   coneP = coneSet(coneP,'spatial density',[0.1 0.5 0.2 0.1]);
%   sensor = sensorCreate('human', coneP);
%   sensorConePlot(sensor)
%
% Copyright ImagEval Consultants, LLC, 2005

if notDefined('sensorName'), sensorName = 'human'; end
if notDefined('pixel'), pixel = pixelCreate('default'); end  % Backward compatibility
    
sensor.name = [];
sensor.type = 'sensor';

% If pixel is really a pixel, then we use its parameters to create a cone
% structure
switch pixel.type
    case 'pixel'
        % Backward compatibility with ISET case
        sensor = sensorSet(sensor,'pixel',pixel);
        sensor = sensorSet(sensor,'spectrum',pixelGet(pixel,'spectrum'));
        sensor = sensorSet(sensor,'size',sensorFormats('qqcif'));
    case 'cone'
        % We sent in a cone, not a pixel
        coneP = pixel;
    otherwise
        error('Bad pixel type, must be pixel or cone')
end

sensor = sensorSet(sensor,'data',[]);
sensor = sensorSet(sensor,'sigmagainfpn',0);    % [V/A]  This is the slope of the transduction function
sensor = sensorSet(sensor,'sigmaoffsetfpn',0);  % V      This is the offset from 0 volts after reset

sensor = sensorSet(sensor,'analogGain',1);
sensor = sensorSet(sensor,'analogOffset',0);
sensor = sensorSet(sensor,'offsetFPNimage',[]);
sensor = sensorSet(sensor,'gainFPNimage',[]);
sensor = sensorSet(sensor,'gainFPNimage',[]);
sensor = sensorSet(sensor,'quantization','analog');

sensorName = ieParamFormat(sensorName);
switch sensorName
    case {'default', 'human'}
        % s = sensorCreate('human', [coneP], [retinalPos], [whichEye]);
        % retinalPos should be 1x2 vector containing eccentricity (deg) and
        % polar angle (deg)
        if notDefined('coneP'), coneP = coneCreate; end
        if ~isempty(varargin)
            retPos = varargin{1};
            if isscalar(retPos), retPos = [retPos, 0]; end
        else
            retPos = [0, 0];
        end
        if length(varargin)>1, whichEye = varargin{2};
        else whichEye = []; end
        
        eccMM = 2*tand(retPos(1)/2) * 17; % assuming focal length of 17 mm
        
        % Assign key fields
        wave = coneGet(coneP,'wave');
        hPixel = pixelCreate('human', wave);
        
        % Adjust pixel gap by retinal position
        coneD = coneDensity(eccMM, retPos(2), whichEye);
        coneSz = sqrt(1/coneD) * 1e-3; % avg cone size with gap in meters
        
        % Adjust pixel gap size according to retinal position
        wGap = coneSz - pixelGet(hPixel, 'width');
        hGap = coneSz - pixelGet(hPixel, 'height');
        assert(wGap>=0 && hGap>=0, 'gap should be non-negative');
        
        hPixel = pixelSet(hPixel, 'width gap', wGap);
        hPixel = pixelSet(hPixel, 'height gap', hGap);
        
        % Add the default human pixel to the sensor
        sensor.spectrum.wave = wave;
        sensor = sensorSet(sensor, 'pixel', hPixel);
        sensor = sensorSet(sensor, 'size', [72 88]);
        
        % Add the default lens structure
        lens = lensCreate([], wave);
        sensor = sensorSet(sensor, 'human lens', lens);
        
        % Add the default macular structure
        macular = macularCreate(macularDensity(retPos(1)), wave);
        sensor = sensorSet(sensor, 'human macular', macular);
        
        % Build up a human cone mosaic.
        sensor = sensorCreateConeMosaic(sensor, coneP);
    case {'color', 'bayer', 'bayer(grbg)', 'bayer-grbg', 'bayergrbg'}
        filterOrder = [2,1;3,2];
        filterFile = 'RGB';
        sensor = sensorBayer(sensor,filterOrder,filterFile);
    case {'bayer(rggb)','bayer-rggb'}
        filterOrder = [1 2 ; 2 3];
        filterFile = 'RGB';
        sensor = sensorBayer(sensor,filterOrder,filterFile);
    case {'bayer(bggr)','bayer-bggr'}
        filterOrder = [3 2 ; 2 1];
        filterFile = 'RGB';
        sensor = sensorBayer(sensor,filterOrder,filterFile);
    case {'bayer(gbrg)','bayer-gbrg'}
        filterOrder = [2 3 ; 1 2];
        filterFile = 'RGB';
        sensor = sensorBayer(sensor,filterOrder,filterFile);

    case 'monochrome'
        % sensorCreate('monochrome')
        filterFile = 'Monochrome';
        sensor = sensorMonochrome(sensor,filterFile);
    case 'monochromearray'
        % nSensors = 5; pixel = [];
        % sensorA = sensorCreate('monochrome array',pixel,nSensors);
        %
        % Builds an array of monochrome sensors, each corresponding to the
        % default monochrome.  The array of sensors is used for
        % calculations that avoid demosaicking.
        if isempty(varargin), N = 3;
        else N = varargin{1};
        end
        
        sensorA(N) = sensorCreate('monochrome');
        for ii=1:(N-1), sensorA(ii) = sensorA(N); end
        sensor = sensorA;
        
        return;
        
    case 'interleaved'
        % Create an interleaved sensor with one transparent and 3 color
        % filters.
        filterFile = 'interleavedRGBW.mat';
        filterPattern = [1 2; 3 4];
        sensor = sensorInterleaved(sensor,filterPattern,filterFile);
    
        
    otherwise
        error('Unknown sensor type');
end

% Set the exposure time - this needs a CFA to be established to account for
% CFA exposure mode.
if sensorCheckHuman(sensor)
    sensor = sensorSet(sensor, 'exp time', 0.05); % 50 ms
    sensor = sensorSet(sensor, 'noise flag', 1); % photons noise only
    sensor = sensorSet(sensor, 'time interval', 0.001);
else
    sensor = sensorSet(sensor,'integrationTime',0);
    sensor = sensorSet(sensor,'autoexposure',1);
    sensor = sensorSet(sensor,'noise flag',2); % all noise
end

sensor = sensorSet(sensor,'CDS',0);

% Put in a default infrared filter.  All ones.
sensor = sensorSet(sensor,'irfilter',ones(sensorGet(sensor,'nwave'),1));

% Place holder for Macbeth color checker positions
sensor = sensorSet(sensor,'mccRectHandles',[]);

return

%-----------------------------
function sensor = sensorBayer(sensor,filterPattern,filterFile)
%
%   Create a default image sensor array structure.

sensor = sensorSet(sensor,'name',sprintf('bayer-%.0f',vcCountObjects('sensor')));
sensor = sensorSet(sensor,'cfa pattern',filterPattern);

% Read in a default set of filter spectra
[filterSpectra,filterNames] = sensorReadColorFilters(sensor,filterFile);
sensor = sensorSet(sensor,'filterspectra',filterSpectra);
sensor = sensorSet(sensor,'filternames',filterNames);

return;

%-----------------------------

%----------------------
% function sensor = sensorMouse(sensor, filterFile)
%
% This isn't right.  The content below should be moved into
% sensorCreateConeMosaic and edited to be made right there.

% error('Not yet implemented');
%
%    sensor = sensorSet(sensor,'name',sprintf('mouse-%.0f',vcCountObjects('sensor')));
%    sensor = sensorSet(sensor,'cfaPattern','mousePattern');
%
%    % try to get the current wavelengths from the scene or the oi.
%    % the mouse sees at different wavelengths than the human : we use
%    % 325-635 usually.
%    scene = vcGetObject('scene');
%    if isempty(scene)
%        getOi = 1;
%    else
%        spect = scene.spectrum.wave;
%        if isempty(spect),  getOi = 1;
%        else
%            mouseWave = spect;
%            getOi = 0;
%        end
%    end
%    if getOi
%       oi = vcGetObject('oi');
%       if isempty(oi), mouseWave = 325:5:635;
%       else spect = oi.optics.spectrum.wave;
%          if isempty(spect),  mouseWave = 325:5:635;
%          else                mouseWave = spect;
%          end
%       end
%    end
%    sensor = sensorSet(sensor,'wave',mouseWave);
%
%    [filterSpectra,filterNames] = sensorReadColorFilters(sensor,filterFile);
%    sensor = sensorSet(sensor,'filterSpectra',filterSpectra);
%    sensor = sensorSet(sensor,'filterNames',filterNames);
%
% return;

%-----------------------------
function sensor = sensorInterleaved(sensor,filterPattern,filterFile)
%
%   Create a default interleaved image sensor array structure.

sensor = sensorSet(sensor,'name',sprintf('interleaved-%.0f',vcCountObjects('sensor')));
sensor = sensorSet(sensor,'cfaPattern',filterPattern);

% Read in a default set of filter spectra
[filterSpectra,filterNames] = sensorReadColorFilters(sensor,filterFile);
sensor = sensorSet(sensor,'filterSpectra',filterSpectra);
sensor = sensorSet(sensor,'filterNames',filterNames);

return;

%-----------------------------
function sensor = sensorCustom(sensor,filterPattern,filterFile)
%
%  Set up a sensor with multiple color filters.
%

sensor = sensorSet(sensor,'name',sprintf('custom-%.0f',vcCountObjects('sensor')));

sensor = sensorSet(sensor,'cfaPattern',filterPattern);

[filterSpectra,filterNames] = sensorReadColorFilters(sensor,filterFile);

% Force the first character of the filter names to be lower case
% This may not be necessary.  But we had a bug once and it is safer to
% force this. - BW
for ii=1:length(filterNames)
    filterNames{ii}(1) = lower(filterNames{ii}(1));
end

sensor = sensorSet(sensor,'filterSpectra',filterSpectra);
sensor = sensorSet(sensor,'filterNames',filterNames);

return;

%-----------------------------
function sensor = sensorMonochrome(sensor,filterFile)
%
%   Create a default monochrome image sensor array structure.
%

sensor = sensorSet(sensor,'name',sprintf('monochrome-%.0f', vcCountObjects('sensor')));

[filterSpectra,filterNames] = sensorReadColorFilters(sensor,filterFile);
sensor = sensorSet(sensor,'filterSpectra',filterSpectra);
sensor = sensorSet(sensor,'filterNames',filterNames);

sensor = sensorSet(sensor,'cfaPattern',1);      % 'bayer','monochromatic','triangle'

return;


