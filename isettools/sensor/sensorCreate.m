function [sensor, params] = sensorCreate(sensorName,pixel,varargin)
%Create an image sensor array structure
%
%   [sensor,params] = sensorCreate(sensorName,[pixel],varargin)
%
% The sensor array uses a pixel definition that can be specified in the
% parameter PIXEL. If this is not passed in, a default PIXEL is created and
% returned.
%
% Several type of image sensors can be created, including multispectral and
% a model of the human cone mosaic.
%
%  Bayer RGB combinations
%      {'bayer-grbg'}
%      {'bayer-rggb'}
%      {'bayer-bggr'}
%      {'bayer-gbrg'}
%
%  Bayer CMY combinations
%      {'bayer (ycmy)'}
%      {'bayer (cyym)'}
%
% Other types
%      {'monochrome'}
%      {'monochrome array'} - sensorCreate('monochrome array',N);
%
% Multiple channel sensors can be created
%      {'grbc'}        - green, red, blue, cyan
%      {'interleaved'} - One transparent channel and 3 RGB.  Same as RGBC
%                        or RGBW
%      {'fourcolor'}
%      {'custom'}
%
% Human cone mosaic
%      {'human'} - Uses Stockman Quanta LMS cones
%                  Default params:  
%                     params.XXX
%
% See also: sensorReadColorFilters, sensorCreateIdeal
%
% Examples
%  sensor = sensorCreate;
%  sensor = sensorCreate('default');
%
%  sensor = sensorCreate('bayer (ycmy)');
%  sensor = sensorCreate('bayer (rggb)');
%  sensor = sensorCreate('Monochrome');
%
%  pSize  = 3e-6;
%  pixel  = [];
%  sensorType = 'rgb';
%  sensor = sensorCreate('ideal',pixel,pSize,sensorType);
%  sensor = sensorCreate('ideal',pixel,pSize,'human','bayer');
%
%  cone   = pixelCreate('human cone'); 
%  sensor = sensorCreate('Monochrome',cone);
%  sensor = sensorCreate('human');
%
%  filterOrder = [1 2 3; 4 5 2; 3 1 4];
%  wave = 400:2:700;
%  filterFile = fullfile(isetRootPath,'data','sensor','colorfilters','sixChannel.mat');
%  pixel = pixelCreate('default',wave);
%  sensorSize = [256 256];
%  sensor = sensorCreate('custom',pixel,filterOrder,filterFile,sensorSize,wave)
%
%  params.sz = [128,192];
%  params.rgbDensities = [0.1 .6 .2 .1]; % Empty (missing cone), L, M, S
%  params.coneAperture = [3 3]*1e-6;     % In meters
%  pixel = [];
%  [sensor, params] = sensorCreate('human',pixel,params);
%  sensorConePlot(sensor)
%
% Copyright ImagEval Consultants, LLC, 2005

if notDefined('sensorName'), sensorName = 'default'; end

sensor.name = [];
sensor.type = 'sensor';

% Make sure a pixel is defined.
if notDefined('pixel')
    pixel  = pixelCreate('default');
    sensor = sensorSet(sensor,'pixel',pixel);
    sensor = sensorSet(sensor,'size',sensorFormats('qqcif'));
else
    sensor = sensorSet(sensor,'pixel',pixel);
end

% The sensor should always inherit the spectrum of the pixel.  Probably
% there should only be one spectrum here, not one for pixel and sensor.
sensor = sensorSet(sensor,'spectrum',pixelGet(pixel,'spectrum'));

sensor = sensorSet(sensor,'data',[]);

sensor = sensorSet(sensor,'sigmagainfpn',0);    % [V/A]  This is the slope of the transduction function
sensor = sensorSet(sensor,'sigmaoffsetfpn',0);  % V      This is the offset from 0 volts after reset

% I wonder if the default spectrum should be hyperspectral, or perhaps it
% should be inherited from the currently selected optical image?
% sensor = initDefaultSpectrum(sensor,'hyperspectral');

sensor = sensorSet(sensor,'analogGain',1);
sensor = sensorSet(sensor,'analogOffset',0);
sensor = sensorSet(sensor,'offsetFPNimage',[]);
sensor = sensorSet(sensor,'gainFPNimage',[]);
sensor = sensorSet(sensor,'gainFPNimage',[]);
sensor = sensorSet(sensor,'quantization','analog');

sensorName = ieParamFormat(sensorName);
switch sensorName
    case {'default','color','bayer','bayer(grbg)','bayer-grbg','bayergrbg'}
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
    case {'bayer(ycmy)','bayer-ycmy'}
        filterFile = 'cym';
        filterOrder = [2,1; 3,2];
        sensor = sensorBayer(sensor,filterOrder,filterFile);
    case {'bayer(cyym)','bayer-cyym'}
        filterFile = 'cym';
        filterOrder = [1 2 ; 2 3];
        sensor = sensorBayer(sensor,filterOrder,filterFile);
    case {'ideal'}
        % sensorCreate('ideal',[],pSize,sensorType,cPattern);
        %
        % sensorType = 'human'  % 'rgb','monochrome'
        % cPattern = 'bayer'    % any sensorCreate option
        % sensorCreate('ideal',[],'human','bayer');
        error('sensorCreate(''ideal'') is deprecated. Set noiseflag to 0');

    case {'custom'}      % Often used for multiple channel
        % sensorCreate('custom',pixel,filterPattern,filterFile,wave);
        if length(varargin) >= 1, filterPattern = varargin{1};
        else  % Must read it here
        end
        if length(varargin) >= 2, filterFile = varargin{2};
        else % Should read it here, NYI
            error('No filter file specified')
        end
        if length(varargin) <= 3 || isempty(varargin{3})
             sensorSize = size(filterPattern);
        else sensorSize = varargin{3};
        end
        if length(varargin) == 4, wave = varargin{4}; 
        else wave = 400:10:700;
        end
        sensor = sensorSet(sensor,'wave',wave);
        sensor = sensorCustom(sensor,filterPattern,filterFile);
        sensor = sensorSet(sensor,'size',sensorSize);
    case {'fourcolor'}  % Often used for multiple channel
        % sensorCreate('custom',pixel,filterPattern,filterFile);
        if length(varargin) >= 1, filterPattern = varargin{1};
        else  % Must read it here
        end
        if length(varargin) >= 2, filterFile = varargin{2};
        else % Should read it here, NYI
            error('No filter file specified')
        end
        sensor = sensorCustom(sensor,filterPattern,filterFile);

    case 'monochrome'
        filterFile = 'Monochrome';
        sensor = sensorMonochrome(sensor,filterFile);
    case 'monochromearray'
        % sensorA = sensorCreate('monochrome array',[],5);
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
    case 'human'
        % sensor = sensorCreate('human',pixel,params);
        % Uses StockmanQuanta
        % See example in header.
        %
        if length(varargin) >= 1, params = varargin{1};
        else params = [];
        end

        % Assign key fields
        if isfield(params,'wave'), wave = params.wave;
        else wave = 400:10:700;
        end
        
        % Add the default human pixel with StockmanQuanta filters.
        sensor = sensorSet(sensor,'wave',wave);
        %         sensor = sensorSet(sensor,'time interval', tInteval);
        sensor = sensorSet(sensor,'pixel',pixelCreate('human',wave));
        
        % Add the default lens structure
        lens = lensCreate([], wave);
        sensor = sensorSet(sensor, 'human lens', lens);
        
        % Add the default macular structure
        macular = macularCreate([], wave);
        sensor = sensorSet(sensor, 'human macular', macular);
             
        % Build up a human cone mosaic.
        sensor = sensorCreateConeMosaic(sensor, params);
        
        % We don't want the pixel to saturate
        pixel  = sensorGet(sensor, 'pixel');
        pixel  = pixelSet(pixel, 'voltage swing', 1);  % 1 volt
        sensor = sensorSet(sensor, 'pixel', pixel);
        
        % There are no filter spectra in the human case.  We calculate the
        % spectral qe from the cones, macular, and lens data
        
    case 'mouse'
        error('NYI: mouse needs to be fixed with sensorCreateConeMosaic');
        %filterFile = 'mouseColorFilters.mat';
        %sensor = sensorMouse(sensor, filterFile);
        %sensor = sensorSet(sensor, 'pixel', pixelCreate('mouse'));
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

return;

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


