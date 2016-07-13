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
%         'couplingFilter',...  - coupling filters time course
%         'couplingMatrix',...  - weights on coupling filters for other cells
%         'generatorFunction',..- the nonlinear function
%         'linearResponse',...  - linear response of all cells
%         'nlResponse',...      - nonlinear response fgenerator(linear) of all cells
%         'numberTrials',...    - number of trials for spike response
%         'spikeResponse',...   - average waveform over N trials including
%                                   post-spike and coupling filter effects
%         'rasterResponse',...  - spike rasters of all cells from N trials
%         'psthResponse'...     - peristimulus time histogram responses of all cells 
% 
% 
% Examples:
%   rgc1.mosaic{1} = mosaicSet(rgc1.mosaic{1}, 'cellType', 'onParasol')
%   rgc1.mosaic{1} = mosaicSet(rgc1.mosaic{1}, 'psthResponse', psth)
% 
% 9/2015 JRG 


% Check for the number of arguments and create parser object.
% Parse key-value pairs.
% 

% % % We could do set using the superclass method
% obj = mosaicSet@rgcMosaic(obj, varargin{:});

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
    'tCenterAll',...
    'tSurroundAll',...
    'responseLinear',...
    'generatorFunction',...
    'nlResponse',...
    'numberTrials',...
    'spikeResponse',...  
    'postSpikeFilter',...
    'couplingFilter',...
    'couplingMatrix',...
    'responseRaster',...
    'responsePsth',...
    'responseVoltage',...
    'responseSpikes'...
    };
p.addRequired('what',@(x) any(validatestring(ieParamFormat(x),allowableFieldsToSet)));
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
switch ieParamFormat(params.what)
    
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
    case{'tcenterall'}
        nCells = size(obj.sRFcenter);
        tCenterNew = cell(nCells(1),nCells(2));
        for ii = 1:nCells(1)
            for jj = 1:nCells(2)
                tCenterNew{ii,jj} = params.value;
            end
        end
        obj.tCenter = tCenterNew;
        
    case{'tsurroundall'}
        nCells = size(obj.sRFsurround);
        tSurroundNew = cell(nCells(1),nCells(2));
        for ii = 1:nCells(1)
            for jj = 1:nCells(2)
                tSurroundNew{ii,jj} = params.value;
            end
        end
        obj.tSurround = tSurroundNew;
    case{'responselinear'}
        obj.responseLinear = params.value;
    case{'generatorfunction'}
        obj.generatorFunction = params.value;        
    case{'nlresponse'}        
        obj.nlResponse = params.value;
    case{'numbertrials'}
        obj.numberTrials = params.value;
    case{'spikeresponse'}
        obj.spikeResponse = params.value;                
    case{'postspikefilter'}
        obj.postSpikeFilter = params.value;
    case{'couplingfilter'}
        obj.couplingFilter = params.value;
    case{'couplingmatrix'}
        obj.couplingMatrix = params.value;
    case{'responseraster'}
        obj.responseRaster = params.value;
    case{'responsepsth'}
        obj.responsePsth = params.value;
    case{'responsevoltage'}
        obj.responseVoltage = params.value;
    case{'responsespikes'}
%         obj.spikeResponse = params.value; 
        nT = size(obj.responseSpikes,3);        
        [sz1,sz2,nTrials,nType] = size(params.value);
%         obj.spikeResponse{1:sz1,1:sz2,nT,1:nType} = params.value;
        if nT == 1 & isempty(obj.responseSpikes); nT = 0; end; 
        
        if nTrials == 1
            for xc = 1:sz1
                for yc = 1:sz2
                    for nTypeI = 1:nType
                        obj.responseSpikes{xc,yc,nT+1,nTypeI} = params.value{xc,yc,1,nTypeI};
                    end
                end
            end
        else
            for xc = 1:sz1
                for yc = 1:sz2
                    for iTrial = 1:nTrials
                        for nTypeI = 1:nType
                            obj.responseSpikes{xc,yc,iTrial,nTypeI} = params.value{xc,yc,iTrial,nTypeI};
                        end
                    end
                end
            end
        end
end

