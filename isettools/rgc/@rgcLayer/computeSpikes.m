function [rgcL, nTrialsSpikeResponse] = computeSpikes(rgcL, varargin)
% IRCOMPUTESPIKES - Generate RGC spikes from the linear mosaic response
%
%   ir = @rgcLayer.computeSpikes(ir, varargin)
%
% Calculate the spikes from the RGC linear response for rgcGLM and rgcLNP
% mosaics.  
%
% Inputs required
%   rgcL:  rgcLayer object
%
% Inputs optiona;
%  nTrialsLinearResponse: Linear inputs from multiple trials
%
% Outputs:
%  The computed spikes are attached to the rgc mosaic
%
% Example:
%   See s_vaRGC.m in WL/WLVernierAcuity
%
% BW (c) isetbio team, 2015

%% Parse
p = inputParser;
p.addRequired('ir',@(x)(isa(rgcL,'rgcLayer')));  % Inner retina object
p.addParameter('coupling',true,@islogical);

p.addParameter('nTrialsLinearResponse',  [], @(x) isnumeric(x)||iscell(x));

p.parse(rgcL,varargin{:});
rgcL = p.Results.ir;
coupling = p.Results.coupling;

nTrialsLinearResponse = p.Results.nTrialsLinearResponse;
%% Required for Pillow code

% The refresh rate refers to the subsampling of the linear response. Each
% sample of the linear response is broken into N bins, where N is the
% refresh rate. This is used in simGLM.m and simGLMcpl.m from Pillow.
% BW - This should be fixed.

global RefreshRate
RefreshRate = 10;

% For every IR, this could be a vector in the future
nRepeats = rgcL.get('number trials');

%% Loop on the mosaics in the inner retina

nTrials = 1;
if ~isempty(nTrialsLinearResponse)
    nTrials = size(nTrialsLinearResponse{1},1);
end

for iTrial = 1:nTrials    
    for ii = 1:length(rgcL.mosaic)
        if ~ismember(class(rgcL.mosaic{ii}), {'rgcGLM','rgcLNP'})
            % No spikes computed
        else
            mosaic   = rgcL.mosaic{ii};
            
            if ~isempty(nTrialsLinearResponse)
                responseLinear = squeeze(nTrialsLinearResponse{ii}(iTrial,:,:,:));
            else
                responseLinear = mosaic.get('response linear');
            end
            nSamples = size(responseLinear,3);
            nCells   = mosaic.get('mosaic size');
            
            % This is a vector of times that each cell spiked
            spikeTimes = cell([nCells,nRepeats]);
            
            % Temporal sample of the voltage response
            respVolts  = zeros(nCells(1),nCells(2),RefreshRate*nSamples,nRepeats);
            
            glminput = RGB2XWFormat(responseLinear);
            
            if coupling
                
                % Call the Pillow code to generate spikes for the whole mosaic
                % using the coupled GLM
                glmprs   = setGLMprs(mosaic,'coupling',coupling);
                
                if ieSessionGet('wait bar'), wbar = waitbar(0,sprintf('Trial 1  of %d',nRepeats)); end
                for tt = 1:nRepeats
                    
                    % Run Pillow code
                    if ieSessionGet('wait bar')
                        wbar = waitbar(0,sprintf('Coupled GLM (Trial %d of %d)',tt,nRepeats));
                    end
                    [responseSpikesVec, Vmem] = simGLMcpl(glmprs, glminput');
                    
                    % Put the data in the outputs
                    cellCtr = 0;
                    for xc = 1:nCells(2)
                        for yc = 1:nCells(1)
                            cellCtr = cellCtr+1;
                            % Vector of times when the cell spiked
                            spikeTimes{yc,xc,tt} = responseSpikesVec{1,cellCtr};
                            respVolts(yc,xc,:,tt)  = Vmem(:,cellCtr);
                        end
                    end
                    if ieSessionGet('wait bar')
                        waitbar(tt/nRepeats,wbar,sprintf('Finished trial %d of %d',tt,nRepeats));
                    end
                end
                
                if ieSessionGet('wait bar'), delete(wbar); end
            else  % No coupling, much faster
                
                % Run without the slow coupling component by looping on the
                % simGLM, not simGLMcp
                
                glmprs   = setGLMprs(mosaic,'coupling',coupling);
                
                for tt = 1:nRepeats
                    
                    cellCtr = 0;   % Reset the cell counter
                    
                    for xc=1:nCells(2)
                        for yc=1:nCells(1)
                            
                            cellCtr = cellCtr+1;
                            % Pull out linear response of a cell
                            thisCellInput = glminput(cellCtr,:);
                            
                            % Run Pillow code
                            [responseSpikesVec, Vmem] = simGLM(glmprs, thisCellInput');
                            
                            % Vector of times when the cell spiked
                            spikeTimes{yc,xc,tt} = responseSpikesVec;
                            
                            % Voltages as a function of time
                            respVolts(yc,xc,:,tt)  = Vmem;
                        end
                    end
                end
                
            end
        end
        
        % Set mosaic spike times
        rgcL.mosaic{ii} = mosaicSet(rgcL.mosaic{ii},'responseSpikes', spikeTimes);
        
        % The nonlinear voltage; only set in the GLM model
        if isa(rgcL.mosaic{ii},'rgcGLM')
            rgcL.mosaic{ii} = mosaicSet(rgcL.mosaic{ii},'responseVoltage', respVolts);
        end
        
        if ~isempty(nTrialsLinearResponse)
            if iTrial == 1 % && ii == 1
                nTrialsSpikeResponse{ii} = ...
                    zeros([nTrials,size(rgcL.mosaic{ii}.get('spikes'))]);
            end
            spikesTemp = rgcL.mosaic{ii}.get('spikes');
            nTrialsSpikeResponse{ii}(iTrial,:,:,1:size(spikesTemp,3)) = spikesTemp;
        end
    end
end

% Maybe this should be [], not 0? (BW)
if isempty(nTrialsLinearResponse), nTrialsSpikeResponse = []; end

end




