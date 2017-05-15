function val = mosaicGet(obj, param, varargin)
% Gets a property from an rgcMosaic object.
% 
%   val = mosaicGet(rgc.mosaic, property)
% 
% The mosaicGet function gets a property from the mosaic object if the
% property is defined in the rgcMosaic superclass. The subclasses of
% rgcMosaic have properties not included here, and if the request is not
% for one of the properties common to all subclasses, then the subclass
% mosaicGet is called.
%
% Inputs: 
%   obj    - rgc object
%   param  - parameter string
%   varargin - Not used yet, but will be used for units and other things.
% 
% Outputs: 
%   val - parameter value
% 
% Properties that can be gotten:
%    'cellType'         - type of RGC of which mosaic is composed
%    'rfDiameter'       - 1 stdev diameter in pixels of spatial RF
%    'rfDiaMagnitude'   - magnitude of spatial RF at 1 stdev
%    'cellLocation'     - coordinates of spatial RF center in ??? frame
%    'sRFcenter'        - center spatial RF surfaces
%    'sRFsurround'      - surround spatial RF surfaces
%    'tCenter'          - center temporal impulse response
%    'tSurround'        - surround temopral impulse response
%    'linearResponse'   - linear response of all cells
%    'dt'
%    'mosaic size'      - Row/col for cells
%    'last spike time'
%    'spikes'
%    'psth'
%
% Examples:
%   val = mosaicGet(rgc1.mosaic{1}, 'cellType')
%   val = mosaicGet(rgc1.mosaic{3}, 'linearResponse')
% 
% 9/2015 JRG (c) isetbio team
% 7/2016 JRG updated
 
%% Parse
p = inputParser; 
p.CaseSensitive = false; 
p.FunctionName = mfilename;
p.KeepUnmatched = true;
 
allowFields = {...
        'celltype',...
        'rfdiameter',...
        'rfdiamagnitude',...
        'celllocation',...
        'srfcenter',...
        'srfsurround',...
        'tcenter',...
        'tsurround',...
        'tonicdrive',...
        'responselinear'...
        'responseSpikes',...
        'spikeTimes',...
        'spikesDownSampled',...
        'mosaicsize', ...
        'dt', ...
        'lastspiketime', ...
        'spikes', ...
        'psth'
    };
p.addRequired('param',@(x) any(validatestring(ieParamFormat(x),allowFields)));
p.addParameter('cell',[],@(x) (length(x(:)) == 2));
 
% Parse and put results into structure p.
p.parse(param,varargin{:}); 
param = ieParamFormat(p.Results.param);
cell = p.Results.cell;
 
%%
val = [];  % Default return value
switch ieParamFormat(param)
    
    case{'celltype'}
        % String that stores cell type name
        val = obj.cellType;
        
    case{'rfdiameter'}
        % Spatial RF diameter in micrometers
        val = obj.rfDiameter; 
        
    case{'rfdiamagnitude'}
        % Magnitude of linear spatial RF in units of spikes/sec at which 1
        % SD contours are computed.
        val = obj.rfDiaMagnitude;
        
    case{'celllocation'}
        % Location of RF center in units of zero-centered cone position
        % The center of the RGC mosaic is at [0 0]
        % If ir.mosaic{1}.cellLocation{1,1} = [-40 -40], then the mosaic
        % underlying cone mosaic is about [80x80] depending on RF size 
        val = obj.cellLocation;
        
    case{'srfcenter'}
        % Linear spatial center RF in units of conditional intensity,
        % related by Poisson firing to spikes/sec.
        val = obj.sRFcenter;
        
    case{'srfsurround'}
        % Linear spatial surround RF in units of conditional intensity,
        % related by Poisson firing to spikes/sec.spikes/sec.
        val = obj.sRFsurround;
        
    case{'tcenter'}
        % Linear temporal center impulse response in units of conditional
        % intensity, related by Poisson firing to spikes/sec
        val = obj.tCenter;
        if ~isempty(cell)
            val = val{cell(1),cell(2)};
        end
        
    case{'tsurround'}
        % Linear temporal surround impulse response in units of conditional
        % intensity, related by Poisson firing to spikes/sec
        val = obj.tSurround;
        if ~isempty(cell)
            val = val{cell(1),cell(2)};
        end
    case{'tonicdrive'}
        % The DC components of the linear response in units of conditional
        % intensity.
        val = obj.tonicDrive;
        if ~isempty(cell)
            val = val{cell(1),cell(2)};
        end
    case{'responselinear'}
        % Linear response in units of conditional intensity, related by
        % Poisson firing to spikes/sec
        val = obj.responseLinear;
        
    case {'mosaicsize'}
        % Mosaic size in units of number of RGCs
        val = size(obj.cellLocation);
        
    case {'dt'}
        % The bin subsampling size. In the original Pillow code, was a
        % fraction of the sampling rate of the linear response (usually
        % 1/120 = 0.0083 sec). Now it takes into account the linear
        % sampling rate and is given in units of microseconds.
        val = obj.dt;   % 10 usec
        
    case {'lastspiketime'}
        % The latest time at which a spike occurs in seconds, used in
        % plotting functions.
        if ~isfield(obj,'responseSpikes')
            return;
        end
        
        nCells  = obj.get('mosaic size');
        nTrials = obj.get('numbertrials');
        spikes  = obj.responseSpikes;
 
        val = 0;
        for ii=1:nCells(1)
            for jj = 1:nCells(2)
                for kk = 1:nTrials
                    mx = max(spikes{ii,jj,kk});
                    val = max([val,mx]);
                end
            end
        end
        
    case{'responsespikes','spiketimes'}
        % Get the spike times in an array
        
        if ~isfield(obj,'responseSpikes')
            return;
        end
        nCells  = obj.get('mosaic size');
        nTrials = obj.get('numbertrials');        
        spikes  = obj.responseSpikes;
        val = 0;
        for ii=1:nCells(1)
            for jj = 1:nCells(2)
                for kk = 1:nTrials
                    mx = length(spikes{ii,jj,kk});
                    val = max([val,mx]);
                end
            end
        end
        val = zeros(nCells(1), nCells(2), mx, nTrials);
        
        for ii=1:nCells(1)
            for jj = 1:nCells(2)
                for kk = 1:nTrials
                    val(ii,jj,1:length(spikes{ii,jj,kk}),kk) = spikes{ii,jj,kk};
                end
            end
        end        
        
    case {'spikes'}
        % What is the difference between this and responseSpikes?
        % cellCtr = 0;
        % @JRG - Needs to be updated
        if ~isfield(obj,'responseSpike')
            return;
        end
        
        dt = obj.dt;
        maxTrials = obj.get('number trials');
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
                    end
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
                
            end
        end
        val = spikes; clear spikes spikesCell;
        
    case{'spikesdownsampled'}
        % Get the spikes in an array, but downsampled by dt
        
        spikes = obj.get('spikes');        
        if isempty(spikes), return; end
        
        nCells  = obj.get('mosaic size');  
                
        numDownSampleBlocks = ceil(size(spikes,3)*obj.dt);
        lenBlock = round(1./obj.dt);
        val = zeros(nCells(1), nCells(2), numDownSampleBlocks);
        for ii=1:nCells(1)
            for jj = 1:nCells(2)
                for kk = 1:numDownSampleBlocks-1
                    downSampleIndStart = (kk-1)*lenBlock + 1;
                    downSampleIndEnd   = (kk)*lenBlock;
                    val(ii,jj,kk) = sum(spikes(ii,jj,downSampleIndStart:downSampleIndEnd));
                end
                kk = numDownSampleBlocks;
                downSampleIndStart = (kk-1)*lenBlock + 1;
                % downSampleIndEnd   = (kk)*lenBlock;
                val(ii,jj,kk) = sum(spikes(ii,jj,downSampleIndStart:end));
            end
        end        
        
    case{'psth'}
        % Calculate the PSTH from the response
        spikes = obj.get('spikes');
        if isempty(spikes), disp('No PSTH'); return; end
        nCells = obj.get('mosaic size');

        dt = obj.dt;
        nSamples = round(1/dt);
        sigma    = nSamples/5;
        convolvewin = fspecial('gaussian',[nSamples,1],sigma);
        % convolvewin = convolvewin/sum(convolvewin(:));
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
        
end
 
end


