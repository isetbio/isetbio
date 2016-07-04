function val = mosaicGet(obj, param, varargin)
% rgcMosaicGet: a method of @rgcMosaic that gets rgcMosaic object
% parameters using the input parser structure.
%
%   val = mosaicGet(rgc.mosaic, param, varargin)
%
% Inputs: 
%   rgc object
%   param - to retrieve
%   vararing depends on parameter
%
% Outputs: 
%  val of parameter
%
% Properties:
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
%         'spikes'
%         'last spike time'
%         'psth'
%
% Examples:
%   val = mosaicGet(rgc1.mosaic{1}, 'cell type')
%   val = mosaicGet(rgc1.mosaic{3}, 'psth response')
%
% 9/2015 JRG

p = inputParser; 
p.CaseSensitive = false; p.FunctionName = mfilename;

% Make key properties that can be set required arguments, and require
% values along with key names.
allowableFields = {...
    'celltype',...
    'rfdiameter',...
    'rfdiaMagnitude',...
    'celllocation',...
    'srfcenter',...
    'srfsurround',...
    'tcenter',...
    'tsurround',...
    'tonicdrive',...
    'responselinear',...
    'generatorfunction',...
    'nlresponse',...
    'postspikefilter',...
    'numbertrials',...
    'responsespikes',...
    'responsevoltage',...
    'responsepsth','psth'...
    'mosaicsize', ...
    'spikes',...
    'lastspiketime'
    };
p.addRequired('param',@(x) any(validatestring(ieParamFormat(x),allowableFields)));
p.addParameter('dt',0.01 , @isnumeric);   % 10 usec

% Parse and put results into structure p.
p.parse(param,varargin{:}); 
param = p.Results.param;
dt    = p.Results.dt;

% Set key-value pairs.
switch ieParamFormat(param)
    case{'celltype'}
        val = obj.cellType;
    case {'mosaicsize'}
        val = size(obj.cellLocation);
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
    case{'tonicdrive'}
        val = obj.tonicDrive;
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
        % Calculate the PSTH from the response
        nCells = obj.get('mosaic size');
        spikes = obj.get('spikes','dt',dt);
                
        nSamples = round(1/dt);
        sigma    = nSamples/5;
        convolvewin = fspecial('gaussian',[nSamples,1],sigma);
        convolvewin = convolvewin/max(convolvewin(:));
        psth = zeros(size(spikes));
        
        for xcell = 1:nCells(1)
            for ycell = 1:nCells(2)
                % The PSTH is the convolved spikes.  We adjusted the window
                % to be sensitive to the selection of dt. Also, we wonder
                % whether the sigma is small enough.  The end points are
                % still noticeably above zero.  We also wonder whether we
                % should normalize, we kind of think not.  BW/HJ.
                
                
                % convolvewin = exp(-(1/2)*(2.5*((0:99)-99/2)/(99/2)).^2);
                % convolvewin= convolvewin./max(convolvewin);
                psth(xcell,ycell,:) = conv(squeeze(spikes(xcell,ycell,:)),convolvewin,'same');
            end
        end
        val = psth;
        
    case {'lastspiketime'}
        maxTrials = obj.numberTrials;
        nCells    = obj.get('mosaic size');
        val = 0;
        for ii=1:nCells(1)
            for jj = nCells(2)
                for kk = 1:maxTrials
                    mx = max(obj.responseSpikes{ii,jj,kk});
                    val = max(val,mx);
                end
            end
        end
        
    case {'spikes'}
        % cellCtr = 0;
        maxTrials = obj.numberTrials;
        nCells    = obj.get('mosaic size');
        lastSpike = obj.get('last spike time');
        spikes = zeros(nCells(1),nCells(2),ceil(lastSpike));
        spikesCell = zeros(maxTrials,ceil(lastSpike/dt));
        
        for xcell = 1:nCells(1)
            for ycell = 1:nCells(2)
                clear spikeTimes spikesCell
                % cellCtr = cellCtr+1;
                
                % Convert the time stamps in response spikes to a vector
                % that has 0's and 1's at different times.  The number
                % of times is 1 ms divided by dt, which appears to be 0.01
                % milliseconds.  So, if there are, say, 500 ms in the
                % stimulus there are 50,000 indices in y
                % We need to speed this up and simplify if possible.
                % 
                for trial = 1:maxTrials
                    spikeTimes =  obj.responseSpikes{xcell,ycell,trial};
                    if strcmpi(class(obj),'rgcphys')
                        % For the rgc physiology in EJ's experiments, the
                        % time base is 10 usec, not 1 ms
                        spikeTimes = .01*spikeTimes;
                    end;
                    % Vector on a time base of dt with a 0 or 1 indicating
                    % a spike or not.
                    spikesCell(trial,ceil(spikeTimes./dt)) = 1;
                end                               

                if size(spikesCell,1) > 1
                    % More than one trial, sum across trials
                    spikes(xcell,ycell,1:length(spikesCell)) = sum(spikesCell);
                else
                    spikes(xcell,ycell,1:length(spikesCell)) = spikesCell;
                end
                val = spikes;
            end
        end
        

end

