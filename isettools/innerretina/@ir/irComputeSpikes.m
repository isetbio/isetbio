function ir = irComputeSpikes(ir, varargin)
% Generate spikes for certain rgc objects
%
% We have models that convert from the continuous response to the spiking
% output for several types of rgc mosaics.  This is the gateway routine
% that examines which rgc type we have and invokes the proper continuous to
% spike response for that type of rgc.
%
% Inputs: the inner retina object
%
% Outputs: the spikes responses are attached to the rgc mosaics in the ir
% object.
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

global RefreshRate
RefreshRate = 100;    
% RefreshRate = exptBinsPerStep


%% Loop on the mosaics in the inner retina
for ii = 1:length(ir.mosaic)
    
    switch class(ir.mosaic{ii})
        case {'rgcGLM','rgcPhys','rgcSubunit'}
            % Call the Pillow code to generate spikes for the whole mosaic
            % using the coupled GLM
            clear responseSpikes responseVoltage
            % Modified
            % responseSpikes = computeSpikesGLM(ir.mosaic{ii,1});
            
            % Wrappers for adapting isetbio mosaic properties to Pillow code
            glminput = setGLMinput(ir.mosaic{ii}.responseLinear);
            glmprs   = setGLMprs(ir.mosaic{ii});
            % Run Pillow code
            [responseSpikesVec, Vmem] = simGLMcpl(glmprs, glminput');
            cellCtr = 0;
            nCells = size(ir.mosaic{ii}.responseLinear);
            responseSpikes = cell(nCells(2),nCells(1));
            responseVoltage = cell(nCells(2),nCells(1));
            for xc = 1:nCells(1)
                for yc = 1:nCells(2)
                    cellCtr = cellCtr+1;
                    responseSpikes{yc,xc} = responseSpikesVec{1,cellCtr};
                    responseVoltage{yc,xc} = Vmem(:,cellCtr);
                end
            end

            % Set mosaic property
            ir.mosaic{ii} = mosaicSet(ir.mosaic{ii},'responseSpikes', responseSpikes);
            ir.mosaic{ii} = mosaicSet(ir.mosaic{ii},'responseVoltage', responseVoltage);
            
        case {'rgcLNP'}
            
            clear responseSpikes responseVoltage
            % This is a place holder for linear, nonlinear, poisson spiking
            % model.  The reference and computations will be explained
            % mainly here.
            glminput = setGLMinput(ir.mosaic{ii}.responseLinear);
            
            % Set the post spike filter to enforce the refractory period.
            glmprs = setPSFprs(ir.mosaic{ii});
            %  glmprs.ih=[]; glmprs.iht=[];
            % No post spike filter - break into different subclass?
            % glmprs = setLNPprs(ir.mosaic{ii});
            
            % Run Pillow code
            if strcmp((ir.name),'Macaque inner retina pixium 1')
                [responseSpikesVec, Vmem] = simGLMcpl(glmprs, glminput');
            else
                [responseSpikesVec, Vmem] = simGLMcpl(glmprs, glminput');
            end
            cellCtr = 0;
            
            nCells = size(ir.mosaic{ii}.responseLinear);
            responseSpikes = cell(nCells(2),nCells(1));
            responseVoltage = cell(nCells(2),nCells(1));
            nCells = size(ir.mosaic{ii}.responseLinear);
            for xc = 1:nCells(1)
                for yc = 1:nCells(2)
                    cellCtr = cellCtr+1;
                    responseSpikes{yc,xc} = responseSpikesVec{1,cellCtr};
                    responseVoltage{yc,xc} = Vmem(:,cellCtr);
                end
            end
            
            ir.mosaic{ii} = mosaicSet(ir.mosaic{ii},'responseSpikes', responseSpikes);
            ir.mosaic{ii} = mosaicSet(ir.mosaic{ii},'responseVoltage', responseVoltage);
            
        otherwise
            error('The rgcMosaic object is a model without a spike response; choose LNP or GLM for spikes.');
    end
end


end



