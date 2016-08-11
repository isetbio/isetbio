function ir = irComputeSpikes(ir, varargin)
% Generate spikes from the linear response
%
%   ir = irComputeSpikes(ir, varargin)
%
% Convert from the linear response to the spiking output for several types
% of rgc mosaics.  Works for rgcGLM and rgcLNP models.
%
% Inputs:
%   inner retina object
%
% Outputs:
%  The spikes responses are attached to the rgc mosaics
%
% Example:
%  os = osCreate('identity');
%  ir = irCreate(os,'model','glm');
%  ir.mosaicCreate('mosaicType','on midget');
%  ir.computeContinuous;
%  ir.computeSpike;
%
% JRG (c) isetbio


%% Parse
p = inputParser;
p.addRequired('ir',@(x)(isa(ir,'ir')));  % Inner retina object
p.addParameter('coupling',true,@islogical);

p.parse(ir,varargin{:});
ir = p.Results.ir;
coupling = p.Results.coupling;

%% Required for Pillow code

% The refresh rate refers to the subsampling of the linear response. Each
% sample of the linear response is broken into N bins, where N is the
% refresh rate. This is used in simGLM.m and simGLMcpl.m from Pillow.
global RefreshRate
RefreshRate = 100;

% For every IR, this could be a vector in the future
nTrials = ir.get('number trials');

%% Loop on the mosaics in the inner retina
for ii = 1:length(ir.mosaic)
    if ~ismember(class(ir.mosaic{ii}), {'rgcGLM','rgcLNP'})
    else
        mosaic = ir.mosaic{ii};
        responseLinear = mosaic.get('response linear');
        nSamples = size(responseLinear,3);
        nCells = mosaic.get('mosaic size');
        
        % This is a vector of times that each cell spiked
        spikeTimes = cell([nCells,nTrials]);
        
        % Temporal sample of the voltage response
        respVolts  = zeros(nCells(1),nCells(2),RefreshRate*nSamples,nTrials);
        
        glminput = RGB2XWFormat(responseLinear);
        
        if coupling
            
            % Call the Pillow code to generate spikes for the whole mosaic
            % using the coupled GLM
            glmprs   = setGLMprs(mosaic,'coupling',coupling);
            
            if ieSessionGet('wait bar'), wbar = waitbar(0,sprintf('Trial 1  of %d',nTrials)); end
            for tt = 1:nTrials
                
                % Run Pillow code
                if ieSessionGet('wait bar')
                    wbar = waitbar(0,sprintf('Coupled GLM (Trial %d of %d)',tt,nTrials));
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
                    waitbar(tt/nTrials,wbar,sprintf('Finished trial %d of %d',tt,nTrials));
                end
            end
            if ieSessionGet('wait bar'), delete(wbar); end
        else
            % Run without the slow coupling component by looping on the
            % simGLM, not simGLMcp
            
            
            glmprs   = setGLMprs(mosaic,'coupling',coupling);
            
            for tt = 1:nTrials
                
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
end

% Set mosaic property
ir.mosaic{ii} = mosaicSet(ir.mosaic{ii},'responseSpikes', spikeTimes);

% The nonlinear voltage which is only set in the GLM model
if isa(ir.mosaic{ii},'rgcGLM')
    ir.mosaic{ii} = mosaicSet(ir.mosaic{ii},'responseVoltage', respVolts);
end

end





