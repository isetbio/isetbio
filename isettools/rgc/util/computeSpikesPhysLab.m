function [spikeTimes spikeDrive psthResponse rollcomp] = computeSpikesPhysLab(obj, varargin)
% Converts the nonlinear response of the generator lookup function to a
% probabilistic spiking output. This is called internally from rgcCompute.
%
% Inputs: the rgc object from rgcCompute.
%
% Outputs:
%       spikeTimes: N trials of spikes computed with the coupled GLM model.
%       spikeDrive: N trials of the drive to the Poisson spike generator.
%       psthResponse: the PSTH combining all N trials of spikes.
%       rollcomp: N trials of the transformed spikeDrive.
%
% Example: see @rgcPhys/rgcCompute.m
%  [spikeResponseFull, spikeDrive, psthResponse, rollcomp] = ...
%                     computeSpikesPhysLab(obj.mosaic{cellTypeInd,1});
%
% (c) isetbio
% 09/2015 JRG

%%%%%% FROM J PILLOW
% -------------  Static nonlinearity & spiking -------------------

ihcpl = mosaicGet(obj, 'couplingFilter');

nCells = length(obj.cellLocation);
nCellsTotal = nCells;%nCells(1)*nCells(2);

Vstm = (obj.responseLinear(1,:,:));
slen = size(Vstm,3);
dt = .1; 

rlen = length([.5+dt:dt:slen-1]');

nbinsPerEval = 100;
nlfun = obj.generatorFunction;
RefreshRate = 120.8;

bindur = dt/RefreshRate;

numberTrials = mosaicGet(obj, 'numberTrials');

% spikeTimes is nCells x 1 in order to make it compatible with other RGC
% implementations where the spatial array is m x n.
spikeTimes = cell(nCells,1,numberTrials,2);

cellCtr = 0;

pairspikecomp = cell(nCells,6,numberTrials);
rollcomp = [];

tic
for xcell = 1:nCells
    rng(1);
    
    % Set coupling filters
    for couplingFilterInd=1:6
        cif_cpgain{couplingFilterInd} = exp(zeros(size((obj.postSpikeFilter{xcell}))));
    end      
    
    Vstm = squeeze(obj.responseLinear(xcell,1,:))';
    slen = length(Vstm);
    
    % Set generator function
    nlfun = obj.generatorFunction{xcell,1};
    if isa(nlfun,'function_handle')
        cif0 = nlfun(reshape( repmat(Vstm, 10, 1) , 1 , slen*10)');
    else
        lnmodel = obj.generatorFunction{xcell,1};
        cif0 = predict(lnmodel, reshape( repmat(Vstm, 10, 1) , 1 , slen*10)');
    end
        
    cif_psgain = exp(obj.postSpikeFilter{xcell});   
    
    ps_bins     = length(cif_psgain);
    cp_bins     = length(cif_cpgain{1});
    
    for i_trial = 1 : numberTrials
 
        cif_ps_cp       = cif0;
        binary_simulation = zeros(1,rlen);
        
        pairspike = zeros(slen*10-190,6) ;
        
        for pair=1:6
            pairspike(pairspikecomp{xcell,pair,i_trial},pair) = 1;
        end

        for i = 1 : slen*10-190%params.bins- max(cp_bins, ps_bins);
            roll = rand(1);
            rollcomp(i_trial,i) = exp(-bindur*cif_ps_cp(i));
            if roll >  exp(-bindur*cif_ps_cp(i));
                cif_ps_cp(i+1: i + ps_bins) =  cif_ps_cp(i+1: i + ps_bins) .* (cif_psgain);
                binary_simulation(i)= 1;
            end
            for pair=1:6%fittedGLM.GLMPars.spikefilters.cp.n_couplings
                
                if pairspike(i,pair) > 0              
                    cif_ps_cp(i+1: i + cp_bins) =  cif_ps_cp(i+1: i + cp_bins) .* (cif_cpgain{pair});
                end
            end
        end
        
        binary_simall{xcell,1,i_trial,1} = binary_simulation; % prune extra zeros
        
        spikeTimes{xcell,1,i_trial,1} = find(binary_simulation==1)'; % prune extra zeros
        spikeTimes{xcell,1,i_trial,2} = cif_ps_cp;
        spikeTimes{xcell,1,i_trial,3} = log(cif0);
        spikeDrive{xcell,1,i_trial} = cif_ps_cp; % prune extra zeros
    end

end
toc

%% Generate raster and PSTH
for ce = 1:nCells;
    
    convolvewin2D = fspecial('gaussian',100,20);
    convolvewin = convolvewin2D(51,:)./max(convolvewin2D(51,:));
    % convolvewin = gausswin(120)/(numberTrials*8.3275e-04*sum(gausswin(120)));
    convolvewin = convolvewin2D(51,:)/(numberTrials*8.3275e-04*sum(convolvewin2D(51,:)));
    
    for trind = 1:numberTrials
        y(trind,:)= binary_simall{ce,1,trind,1};
    end
    
    PSTH_rec=conv(sum(y),convolvewin,'same');
    
    psthResponse{ce} = PSTH_rec;
end


end


