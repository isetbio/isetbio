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

%% Required for Pillow code

% To be eliminated

% fps = 1/125;        % frame rate of 125 fps; usually 121 in lab but 125 for integer steps
% normalRR = fps;
% 
% exptRR = irGet(ir,'timing');
% exptBinsPerStep = round(normalRR/exptRR);

% @JRG - To comment
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

        spikeTimes = cell([nCells,nTrials]);
        respVolts  = zeros(nCells(1),nCells(2),RefreshRate*nSamples,nTrials);
        
        % Call the Pillow code to generate spikes for the whole mosaic
        % using the coupled GLM
        glminput = RGB2XWFormat(responseLinear);
        glmprs   = setGLMprs(mosaic);
        
        for tt = 1:nTrials
    
            % Run Pillow code
            [responseSpikesVec, Vmem] = simGLMcpl(glmprs, glminput');
            
            % Put the data in the outputs
            cellCtr = 0;
            for xc = 1:nCells(1)
                for yc = 1:nCells(2)
                    cellCtr = cellCtr+1;
                    spikeTimes{yc,xc,:,tt} = responseSpikesVec{1,cellCtr};
                    respVolts(yc,xc,:,tt)  = Vmem(:,cellCtr);
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
end

end



