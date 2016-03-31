function val = mosaicGet(obj, varargin)
% rgcMosaicGet: a method of @rgcMosaic that gets rgcMosaic object 
% parameters using the input parser structure.
% 
%       val = mosaicGet(rgc.mosaic, property)
% 
% Inputs: rgc object, property to be gotten
% 
% Outputs: val of property
% 
% Properties that can be gotten:
%         'cellType',...        - type of RGC of which mosaic is composed
%         'rfDiameter',...      - 1 stdev diameter in pixels of spatial RF
%         'rfDiaMagnitude',...  - magnitude of spatial RF at 1 stdev
%         'cellLocation',...    - coordinates of center of spatial RF
%         'sRFcenter',...       - center spatial RF surfaces
%         'sRFsurround',...     - surround spatial RF surfaces
%         'tCenter',...         - center temporal impulse response
%         'tSurround',...       - surround temopral impulse response
%         'postSpikeFilter',... - post-spike filter time course
%         'couplingFilter',...  - coupling filters time course
%         'couplingMatrix',...  - weights on coupling filters for other cells
%         'generatorFunction',..- the nonlinear function
%         'responseLinear',...  - linear response of all cells
%         'nlResponse',...      - nonlinear response fgenerator(linear) of all cells
%         'numberTrials',...    - number of trials for spike response
%         'responseSpikes',...   - average waveform over N trials including
%                                   post-spike and coupling filter effects
%         'responseVoltage',... - the voltage waveform used for coupling 
% 
% 
% Examples:
%   val = mosaicGet(rgc1.mosaic{1}, 'cellType')
%   val = mosaicGet(rgc1.mosaic{3}, 'psthResponse')
% 
% 9/2015 JRG 

% Check key names with a case-insensitive string, errors in this code are
% attributed to this function and not the parser object.
error(nargchk(0, Inf, nargin));
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
    'responseLinear',...
    'generatorFunction',...
    'numberTrials',...
    'responseSpikes',...   
    'responseVoltage',...
    'couplingFilter',...
    'couplingMatrix',...   
    'postSpikeFilter',...
    'numberSubunits'...
    };
p.addRequired('what',@(x) any(validatestring(x,allowableFieldsToSet)));

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
        val = obj.cellType;
    case{'rfdiameter'}
        val = obj.rfDiameter;
    case{'rfdiamagnitude'}
        val = obj.rfDiaMagnitude;
    case{'celllocation'}
        val = obj.cellLocation;
    case{'srfcenter'}
        val = obj.sRFcenter;
    case{'srfsurround'}
        val = obj.sRFsurround;
    case{'tcenter'}
        val = obj.tCenter;
    case{'tsurround'}
        val = obj.tSurround;
    case{'responselinear'}
        val = obj.responseLinear;
    case{'generatorfunction'}
        val = obj.generatorFunction;
    case{'responsespikes'}
        val = obj.responseSpikes;
    case{'responsevoltage'}
        val = obj.responseVoltage;        
    case{'numbertrials'}
%         val = obj.numberTrials;
        val = size(obj.responseSpikes,3);
    case{'postspikefilter'}
        val = obj.postSpikeFilter;
    case{'couplingfilter'}
        val = obj.couplingFilter;
    case{'couplingmatrix'}
        val = obj.couplingMatrix;
    case{'numbersubunits'}
        val = obj.numberSubunits;
    case{'rasterresponse'}
        val = obj.rasterResponse;
    case{'psthresponse'}
        val = obj.psthResponse;
end

