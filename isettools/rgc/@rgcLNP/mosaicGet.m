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
%         'generatorFunction',..- the nonlinear function
%         'responseLinear',...  - linear response of all cells
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
    'responseLinear',...
    'generatorFunction',...
    'nlResponse',...
    'postSpikeFilter',...
    'numberTrials',...
    'responseSpikes',...    
    'responseVoltage',...    
    'responsePsth','psth'...
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
    case{'postspikefilter'}
        val = obj.postSpikeFilter;
    case{'numbertrials'}
        % val = obj.numberTrials;       
        val = size(obj.responseSpikes,3);
    case{'responsespikes'}
        val = obj.responseSpikes;
    case{'responsevoltage'}
        val = obj.responseVoltage;
    case{'couplingfilter'}
        val = obj.couplingFilter;
    case{'couplingmatrix'}
        val = obj.couplingMatrix;        
    case{'responsepsth', 'psth'}
 
% %         if ~isempty(obj.responsePsth)
% %             val = obj.responsePsth;
% %         else
%             nCells = size(obj.responseSpikes,1)*size(obj.responseSpikes,2);
%             numberTrials = obj.numberTrials;
%             for ce = 1:nCells
%                 convolvewin2D = fspecial('gaussian',100,20);
%                 convolvewin = convolvewin2D(51,:)./max(convolvewin2D(51,:));
% %                 convolvewin=gausswin(100);
% %                 y = zeros(numberTrials,100000);
%                 clear y
%                 for trind = 1:numberTrials
%                     y(trind,ceil(100*obj.responseSpikes{ce,1,trind}))=1;
%                 end
%                 PSTH_rec=conv(sum(y),convolvewin,'same');               
%                 psthResponse{ce} = PSTH_rec;
%             end
%             clear y
%             val = psthResponse;
% %         end

        cellCtr=0; dt = .01;
        maxTrials = obj.numberTrials;
        nCells = size(obj.responseSpikes);
        yout = [];
            for xcell = 1:nCells(1)
                for ycell = 1:nCells(2)
                    clear yind y
                    cellCtr = cellCtr+1;
                    
                    for trial = 1:maxTrials
                        
                        yind =  obj.responseSpikes{xcell,ycell,trial,1};
                        if strcmpi(class(obj),'rgcphys'); yind = .01*yind; end;
                        y(trial,ceil(yind./dt))=1;
                    end
    
                    [jv,iv] = ind2sub([nCells(1),nCells(2)],cellCtr);
                    cellCtr2 = sub2ind([nCells(2),nCells(1)],iv,jv);
                    
                    %                     subplot(nCells(2),nCells(1),cellCtr);
                    
                    %
                    convolvewin = exp(-(1/2)*(2.5*((0:99)-99/2)/(99/2)).^2);                                       
                    convolvewin= convolvewin./max(convolvewin);
                    yout(cellCtr,1:length(y)) = y;%sum(y);
                    PSTH_out{xcell,ycell}=conv(sum(y),convolvewin,'same');
%                     plot(.1*bindur:.1*bindur:.1*bindur*length(PSTH_rec),PSTH_rec);
                end
            end
            
            val.psth = PSTH_out;
            val.spikes = yout;
            
        end
end

