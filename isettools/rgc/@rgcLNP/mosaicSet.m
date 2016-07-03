function obj = mosaicSet(obj, varargin)
% rgcMosaicSet: a method of @rgcMosaic that sets rgcMosaic object 
% parameters using the input parser structure.
% 
%       rgc.mosaic = mosaicSet(rgc.mosaic, property, value)
%  
% Inputs: rgc object, key-value pair of property and value to which it is
%   being set.
% 
% Outputs: obj with property set appropriately.% 
% 
% Properties that can be set:
%         'cellType',...        - type of RGC of which mosaic is composed
%         'rfDiameter',...      - 1 stdev diameter in pixels of spatial RF
%         'rfDiaMagnitude',...  - magnitude of spatial RF at 1 stdev
%         'cellLocation',...    - coordinates of center of spatial RF
%         'sRFcenter',...       - center spatial RF surfaces
%         'sRFsurround',...     - surround spatial RF surfaces
%         'tCenter',...         - center temporal impulse response
%         'tSurround',...       - surround temopral impulse response
%         'postSpikeFilter',... - post-spike filter time course
%         'generatorFunction',..- the nonlinear function
%         'responseLinear',...  - linear response of all cells
%         'numberTrials',...    - number of trials for spike response
%         'responseSpikes',...  - average waveform over N trials including
%                                   post-spike and coupling filter effects
%         'responseVoltage',... - the voltage waveform used for coupling 
% 
% Examples:
%   rgc1.mosaic{1} = mosaicSet(rgc1.mosaic{1}, 'cellType', 'onParasol')
%   rgc1.mosaic{1} = mosaicSet(rgc1.mosaic{1}, 'psthResponse', psth)
% 
% 9/2015 JRG 
% 
% 9/2015 JRG 

% % % We could do set using the superclass method
% obj = mosaicSet@rgcMosaic(obj, varargin{:});

% Check key names with a case-insensitive string, errors in this code are
% attributed to this function and not the parser object.
narginchk(0, Inf);
p = inputParser; p.CaseSensitive = false; p.FunctionName = mfilename;

% Make key properties that can be set required arguments, and require
% values along with key names.
allowableFieldsToSet = {...
    'cellType',...
    'rfDiameter',...
    'rfDiaMagnitude',...
    'cellLocation',...
    'sRFcenter',...
    'sRFsurround',...
    'tCenter',...
    'tSurround',...
    'tonicDrive',...
    'responseLinear',...
    'generatorFunction',...
    'nlResponse',...
    'postSpikeFilter',...
    'numberTrials'...
    'responseSpikes',...  
    'responseVoltage',...
    'rasterResponse',...
    'psthResponse'...
    };
p.addRequired('what',@(x) any(validatestring(x,allowableFieldsToSet)));
p.addRequired('value');

% % Define what units are allowable.
% allowableUnitStrings = {'a', 'ma', 'ua', 'na', 'pa'}; % amps to picoamps
% 
% % Set up key value pairs.
% % Defaults units:
% p.addParameter('units','pa',@(x) any(validatestring(x,allowableUnitStrings)));

% Parse and put results into structure p.
p.parse(varargin{:}); params = p.Results;

% % Old error check on input.
% if ~exist('params','var') || isempty(params)
%     error('Parameter field required.');
% end
% if ~exist('val','var'),   error('Value field required.'); end;

% Set key-value pairs.
switch lower(params.what)
    
    case{'celltype'}
        obj.cellType = params.value;
    case{'rfdiameter'}
        obj.rfDiameter = params.value;
    case{'rfdiamagnitude'}
        obj.rfDiaMagnitude = params.value;
    case{'celllocation'}
        obj.cellLocation = params.value;
    case{'srfcenter'}
        obj.sRFcenter = params.value;
    case{'srfsurround'}
        obj.sRFsurround = params.value;
    case{'tcenter'}
        obj.tCenter = params.value;
    case{'tsurround'}
        obj.tSurround = params.value;
    case{'tonicdrive'}
        for ci1 = 1:size(obj.tonicDrive,1)
            for ci2 = 1:size(obj.tonicDrive,2)
                
                obj.tonicDrive{ci1,ci2} = params.value;
            end
        end
    case{'responselinear'}
        obj.responseLinear = params.value;
    case{'generatorfunction'}
        obj.generatorFunction = params.value;   
    case{'postspikefilter'}
        obj.postSpikeFilter = params.value;
    case{'numbertrials'}
        obj.numberTrials = params.value;
    case{'responsespikes'}
%         obj.spikeResponse = params.value; 
        nT = size(obj.responseSpikes,3);        
        [sz1,sz2,nTrials,nType] = size(params.value);
%         obj.spikeResponse{1:sz1,1:sz2,nT,1:nType} = params.value;
        if nT == 1 & isempty(obj.responseSpikes); nT = 0; end; 
        for xc = 1:sz1
            for yc = 1:sz2
                for nTypeI = 1:nType
                    obj.responseSpikes{xc,yc,nT+1,nTypeI} = params.value{xc,yc,1,nTypeI};
                end
            end
        end      
    case{'responsevoltage'}
        obj.responseVoltage = params.value;
end

