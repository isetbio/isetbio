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
%   val = mosaicGet(rgc1.mosaic{1}, 'cellType')
%   val = mosaicGet(rgc1.mosaic{3}, 'psthResponse')
% 
% 9/2015 JRG 

% Check for the number of arguments and create parser object.
% Parse key-value pairs.
% 

% % % We could do set using the superclass method
% val = mosaicGet@rgcMosaic(obj, varargin{:});

% Check key names with a case-insensitive string, errors in this code are
% attributed to this function and not the parser object.
narginchk(0, Inf);
p = inputParser; p.CaseSensitive = false; p.FunctionName = mfilename;
p.KeepUnmatched = true;
% Make key properties that can be set required arguments, and require
% values along with key names.
allowableFieldsToSet = {...
    'cellType',...
    'rfDiameter',...
    'rfDiaMagnitude',...
    'cellLocation',...
    'mosaicSize',...
    'sRFcenter',...
    'sRFsurround',...
    'tCenter',...
    'tSurround',...
    'tonicDrive',...
    'linearResponse',...
    'responseLinear',...
    'generatorFunction',...
    'nlResponse',...
    'numberTrials',...
    'spikeResponse',... 
    'postSpikeFilter',...
    'couplingFilter',...
    'couplingMatrix',...
    'lastspiketime',....
    'responseRaster',...
    'responsePsth',...
    'responseVoltage',...
    'responseSpikes'...
    };
p.addRequired('what',@(x) any(validatestring(ieParamFormat(x),allowableFieldsToSet)));
p.addParameter('cell',[],@(x) (length(x(:)) == 2));

% % Define what units are allowable.
% allowableUnitStrings = {'a', 'ma', 'ua', 'na', 'pa'}; % amps to picoamps
% 
% % Set up key value pairs.
% % Defaults units:
% p.addParameter('units','pa',@(x) any(validatestring(x,allowableUnitStrings)));

% Parse and put results into structure p.
p.parse(varargin{:}); params = p.Results;
cell = p.Results.cell;

% % Old error check on input.
% if ~exist('params','var') || isempty(params)
%     error('Parameter field required.');
% end
% if ~exist('val','var'),   error('Value field required.'); end;

% Set key-value pairs.
switch ieParamFormat(params.what)
    case{'celltype'}
        val = obj.cellType;
    case{'rfdiameter'}
        val = obj.rfDiameter;
    case{'rfdiamagnitude'}
        val = obj.rfDiaMagnitude;
    case{'celllocation'}
        val = obj.cellLocation;
        
    case {'mosaicsize'}
        % Mosaic size in units of number of RGCs
        val = size(obj.cellLocation);
    case{'srfcenter'}
        val = obj.sRFcenter;
    case{'srfsurround'}
        val = obj.sRFsurround;
    case{'tcenter'}
        val = obj.tCenter;
        if ~isempty(cell)
            val = val{cell(1),cell(2)};
        end
    case{'tsurround'}
        val = obj.tSurround;
        if ~isempty(cell)
            val = val{cell(1),cell(2)};
        end
    case{'tonicdrive'}
        val = obj.tonicDrive;        
        if ~isempty(cell)
            val = val{cell(1),cell(2)};
        end
    case{'linearresponse'}
        val = obj.linearResponse;
    case{'responselinear'}
        val = obj.responseLinear;
    case{'generatorfunction'}
        val = obj.generatorFunction;
    case{'nlresponse'}
        val = obj.nlResponse;
    case{'numbertrials'}
        val = obj.numberTrials;
    case{'spikeresponse'}
        val = obj.spikeResponse;
    case{'postspikefilter'}
        val = obj.postSpikeFilter;
    case{'couplingfilter'}
        val = obj.couplingFilter;
    case{'couplingmatrix'}
        val = obj.couplingMatrix;
        
    case {'lastspiketime'}
        nCells    = obj.get('mosaic size');
        nTrials = obj.get('numbertrials');
        spikes = obj.responseSpikes;

        val = 0;
        for ii=1:nCells(1)
            for jj = 1:nCells(2)
                for kk = 1:nTrials
                    mx = max(spikes{ii,jj,kk});
                    val = max([val,mx]);
                end
            end
        end
    case{'responseraster'}
        val = obj.responseRaster;
    case{'responsepsth'}
        if ~isempty(obj.responsePsth)
            val = obj.responsePsth;
        else
            % convolve = gausswin(smoothing)/(trials*fittedLN.orig_glm.t_bin*sum(gausswin(smoothing)));
            % LN_FR = conv(sum(fittedLN.xval.rasters.ln_sim),convolve,'same');

            nCells = size(obj.responseSpikes,1)*size(obj.responseSpikes,2);
            numberTrials = obj.numberTrials;
            for ce = 1:nCells
%                 convolvewin2D = fspecial('gaussian',100,20);
%                 convolvewin = convolvewin2D(51,:)./max(convolvewin2D(51,:));
%                 convolvewin=gausswin(100);
                convolvewin = gausswin(120)/(numberTrials*8.3275e-04*sum(gausswin(120)));
                clear y
                for trind = 1:numberTrials
                    y(trind,obj.responseSpikes{ce,1,trind})=1;
                end
                PSTH_rec=conv(sum(y),convolvewin,'same');               
                psthResponse{ce} = PSTH_rec;
            end
            val = psthResponse;
        end
        
    case{'responsevoltage'}
        val = obj.responseVoltage;
end

  