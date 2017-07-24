function val = get(obj, param, varargin)
% Gets a property from an rgcMosaic object.
% 
%   val = rgcMosaic.get(rgc.mosaic, property)
% 
% Gets a property from an rgc mosaic object. 
%
% The subclasses of rgcMosaic (e.g., rgcGLM and rgcLNP) have properties not
% included here. Those classes first check their own special cases, and if
% the request is not for one of those properties then this function is
% called for a more general case common to all rgc subclasses.
%
% Inputs: 
%   obj      - rgc object
%   param    - parameter string
%   varargin - This will be used for units and other things.
% 
% Outputs: 
%   val - parameter value
% 
% Use @rgcMosaic.get('help') to see allowable parameters.
%
% Examples:
%   val = @rgcMosaic.get('cellType')
%   val = @rgcMosaic.get('linearResponse')
% 
% 9/2015 JRG, BW (c) isetbio team
 
%% Parse
p = inputParser; 
p.CaseSensitive = false; 
p.FunctionName = mfilename;
p.KeepUnmatched = true;
 
allowFields = {...
        'help',...
        'celltype',...
        'ntrials',...
        'rfdiameter',...
        'rfdiamagnitude',...
        'celllocation',...
        'srfcenter',...
        'srfsurround',...
        'spatialrf', ...
        'tcenter',...
        'tsurround',...
        'tonicdrive',...
        'responselinear'...
        'spiketimes',...
        'spikesdownsampled',...
        'mosaicsize', ...
        'mosaicsamples', ...
        'meterspersample', ...
        'dt', ...
        'lastspiketime', ...
        'spikes', ...
        'psth'
    };
p.addRequired('param',@(x) any(validatestring(ieParamFormat(x),allowFields)));
 
% Parse the arguments into the p.Results field
p.parse(param,varargin{:}); 

% cell = p.Results.cell;
 
val = [];  % Default return value

%% Remove the spaces and force lower case on the param argument
switch ieParamFormat(param)
        case 'help'
        fprintf('\nKnown %s parameters\n--------------\n',class(obj));
        for ii=2:length(allowFields)
            fprintf('\t%s\n',allowFields{ii});
        end
        return;
    case{'celltype'}
        % String that stores cell type name
        val = obj.cellType;
    case {'ntrials'}
        val = obj.parent.nTrials;
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
        
    case {'spatialrf'}
        val = obj.sRFcenter - obj.sRFsurround;
        
    case{'tcenter'}
        % Linear temporal center impulse response in units of conditional
        % intensity, related by Poisson firing to spikes/sec
        val = obj.tCenter;
        %         if ~isempty(cell)
        %             val = val{cell(1),cell(2)};
        %         end
        
    case{'tsurround'}
        % Linear temporal surround impulse response in units of conditional
        % intensity, related by Poisson firing to spikes/sec
        val = obj.tSurround;
        %         if ~isempty(cell)
        %             val = val{cell(1),cell(2)};
        %         end
    case{'tonicdrive'}
        % The DC components of the linear response in units of conditional
        % intensity.
        val = obj.tonicDrive;
        %         if ~isempty(cell)
        %             val = val{cell(1),cell(2)};
        %         end
    case{'responselinear'}
        % Linear response in units of conditional intensity, related by
        % Poisson firing to spikes/sec
        val = obj.responseLinear;
        
    case {'mosaicsize'}
        val = obj.patchSize;
        if isempty(val), val = obj.input.patchSize; end
        
    case {'mosaicsamples'}        
        % Mosaic size in units of number of RGCs
        [nRow, nCol, ~] = size(obj.cellLocation);
        val = [nRow, nCol];
        
    case {'meterspersample'}
        % obj.get('meters per sample')
        
        % Find the range of the cell locations.  This corresponds to the
        % number of samples on the input layer, really.
        r = obj.cellLocation(:,:,1); r = r(:); rRange = range(r);
        c = obj.cellLocation(:,:,1); c = c(:); cRange = range(c);
        
        % Find the patch size.
        patchSize = obj.get('mosaic size');     % Meters
        
        % Divide
        val = patchSize ./ [rRange, cRange] ;
        
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
        
        nCells  = obj.get('mosaic samples');
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
        
    case{'spiketimes'}
        % Get the spike times in an array
        
        if ~isfield(obj,'responseSpikes')
            return;
        end
        nCells  = obj.get('mosaic samples');
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
        % Returns 
        
        dt = obj.dt;  % Bin time in microseconds
        
        maxTrials = obj.get('n trials');
        nCells    = obj.get('mosaic samples');
        lastSpike = obj.get('last spike time');
        spikes = zeros(nCells(1),nCells(2),ceil(lastSpike));
        spikesCell = zeros(maxTrials,ceil(lastSpike/dt));
        
        for xcell = 1:nCells(1)
            for ycell = 1:nCells(2)
                clear spikeTimes spikesCell
                
                % Convert the time stamps in response spikes to a vector
                % that has 0's and 1's at different times.  The number of
                % sample times is 1 ms divided by dt, which appears to be
                % 0.1 (100 usec) or 0.01 milliseconds (10 usec)
                %
                % So, if there are, say, 500 ms in the stimulus there are
                % 50,000 indices in spikeTimes. 
                %
                for trial = 1:maxTrials
                    % These times are every dt*msec step.  Typically dt is
                    % 0.1 or 0.01.
                    spikeTimes =  obj.responseSpikes{xcell,ycell,trial};
                    %   if strcmpi(class(obj),'rgcphys')
                    %       % For the rgc physiology in EJ's experiments, the
                    %       % time base is 10 usec, not 1 ms
                    %      spikeTimes = .01*spikeTimes;
                    %   end
                    
                    % This is a vector of 0s and 1s on a time base of
                    % milliseconds.
                    spikesCell(trial,ceil(spikeTimes./dt)) = 1;
                end
                
                if size(spikesCell,1) > 1
                    % More than one trial, sum across trials
                    disp('Surprised to be rgcMosaic.mosaicGet in this case')
                    spikes(xcell,ycell,1:length(spikesCell)) = sum(spikesCell);
                else
                    spikes(xcell,ycell,1:length(spikesCell)) = spikesCell;
                end
                
            end
        end
        val = spikes; 
        clear spikes spikesCell;
        
    case{'spikesdownsampled'}
        % Get the spikes in an array, but downsampled by dt
        
        spikes = obj.get('spikes');        
        if isempty(spikes), return; end
        
        nCells  = obj.get('mosaic samples');  
                
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
        nCells = obj.get('mosaic samples');

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


