
function [spikeTimes spikeDrive psthResponse rollcomp] = computeSpikesPhysLab(obj, varargin)
% computeSpikes: a util function of the @rgc parent class, this
% converts the nonlinear response of the generator lookup function to a
% probabilistic spiking output.
%
% Inputs:
%
% Outputs:
%
% Example:
%
% (c) isetbio
% 09/2015 JRG

%%%%%% FROM J PILLOW
% -------------  Static nonlinearity & spiking -------------------

% ih = mosaicGet(obj, 'postSpikeFilter');
% rng(1); % move this?
ihcpl = mosaicGet(obj, 'couplingFilter');

spResponseSize = size(obj.nlResponse{1,1}(:,:,1));
nSamples = size(obj.nlResponse{1,1},3);

nCells = length(obj.cellLocation);
nCellsTotal = nCells;%nCells(1)*nCells(2);

Vstm = obj.nlResponse{1,1};
slen = length(Vstm);
dt = .1; % sensorGet(sensor,'integration time');
% rlen = length([.5+dt:dt:slen+.5]');
rlen = length([.5+dt:dt:slen-1]');
% hlen = length(ihcpl{1}{1});
nbinsPerEval = 100;
nlfun = obj.generatorFunction;
RefreshRate = 120;

% bindur = dt/RefreshRate;
bindur = 1/1208;
numberTrials = mosaicGet(obj, 'numberTrials');

% spikeTimes is nCells x 1 in order to make it compatible with other RGC
% implementations where the spatial array is m x n.
spikeTimes = cell(nCells,1,numberTrials,2);

cellCtr = 0;

% load('/Users/james/Documents/matlab/rgc Parameters/pairspike_1.mat','pairspike');
load('/Users/james/Documents/matlab/rgc Parameters/pairspikeall.mat','pairspike');

nlfun = obj.generatorFunction;
tic
for xcell = 1:nCells
    rng(1);
    for couplingFilterInd=1:6
        cif_cpgain{couplingFilterInd} = exp(obj.couplingFilter{xcell}{couplingFilterInd});
    end
    cif_psgain = exp(obj.postSpikeFilter{xcell});
    ps_bins     = length(cif_psgain);
    cp_bins     = length(cif_cpgain{1});
    
    Vstm = vertcat(obj.linearResponse{:,xcell,1});
    slen = length(Vstm);
    % cif0 = nlfun(interp1([0:slen-1]',Vstm',[.5+dt:dt:slen-1]', 'linear'));
    cif0 = nlfun(reshape( repmat(Vstm, 10, 1) , 1 , 6300)');
    % lcif_kx0 = reshape( repmat(lcif_kx_frame, bpf, 1) , 1 , params.bins);
    for i_trial = 1 : numberTrials
%         tic    
        cif_ps_cp       = cif0;
        binary_simulation = zeros(1,rlen);
        for i = 1 : 6010%params.bins- max(cp_bins, ps_bins);
            roll = rand(1);
            rollcomp(i_trial,i) = exp(-bindur*cif_ps_cp(i));
            if roll >  exp(-bindur*cif_ps_cp(i));
                cif_ps_cp(i+1: i + ps_bins) =  cif_ps_cp(i+1: i + ps_bins) .* (cif_psgain);
                binary_simulation(i)= 1;
            end
            for pair=1:6%fittedGLM.GLMPars.spikefilters.cp.n_couplings
                if pairspike{xcell,pair}(i_trial,i)
                    cif_ps_cp(i+1: i + cp_bins) =  cif_ps_cp(i+1: i + cp_bins) .* (cif_cpgain{pair});
                end
            end
        end
%         logical_sim(i_trial,:) = binary_simulation ;
%         drive(i_trial, :) = bindur*cif_ps_cp;
%         toc
        spikeTimes{xcell,1,i_trial,1} = binary_simulation; % prune extra zeros
        spikeTimes{xcell,1,i_trial,2} = cif0;
        spikeTimes{xcell,1,i_trial,3} = log(cif0);
        spikeDrive{xcell,1,i_trial} = cif_ps_cp; % prune extra zeros
    end
    toc
end
toc
ph=1;

%%

figure; 
hold on; for ce = 1:nCells; 
    
%     for tr = 1:numberTrials;
%         %       clear yind y
%         % subplot(6,6,ce);
%         % if ~isempty(spikeTimes{ce,1,tr,1});
%         % subplot(6,7,ce); hold on; plot(spikeTimes{ce,1,tr,1},tr,'ok');axis([0 270 0 10]);end;end;end;
%         subplot(2,1,1);
%         spikeTimesP = find(spikeTimes{ce,1,tr,1} == 1);
%         hold on; line([spikeTimesP',spikeTimesP'].*bindur,[tr tr-1],'color','k');
%         axis([0 5000 0 numberTrials]);
%         
%         % end;
%     end%trials;
%     % end;



% subplot(2,1,2);
% subplot(6,7,ce);
convolvewin=gausswin(100);
for trind = 1:numberTrials
%     yind= spikeTimes{ce,1,trind,1};
%     y(trind,round(yind./dt))=1;

    y(trind,:)= spikeTimes{ce,1,trind,1};
end
% yind=sort(vertcat(spikeTimes{ce,1,1:numberTrials,1}),'ascend');
% y(round(yind./dt))=1;
PSTH_rec=conv(sum(y),convolvewin,'same');
plot(bindur:bindur:bindur*length(PSTH_rec),PSTH_rec);        
axis([0 500 0  max(PSTH_rec)])
[maxv, maxi] = max(PSTH_rec); title(sprintf('maxv = %.1f, maxi = %d',maxv,maxi));
%     end;

psthResponse{ce} = PSTH_rec;
end

maxVal = 100;%max(abs(horzcat(meanVoltage{:})));
axesHandles = get(gcf,'children');
if isnan(maxVal), maxVal = 0.00001; end;
axis(axesHandles,[0 5 0 maxVal])
clear axesHandles;

end


