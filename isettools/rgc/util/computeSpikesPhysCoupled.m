
function spikeTimes = computeSpikesPhysCoupled(obj, varargin)
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

% ihthi = [dt:dt:max(glmprs.iht)]';  % time points for sampling
% ihthi = [dt:dt:slen*dt];
% ihhi = interp1(.001:.01:slen*dt+.01, ih, ihthi, 'linear', 0);
% hlen = length(ihhi);

numberTrials = mosaicGet(obj, 'numberTrials');

% spikeTimes is nCells x 1 in order to make it compatible with other RGC
% implementations where the spatial array is m x n.
spikeTimes = cell(nCells,1,numberTrials,2);

cellCtr = 0;

for xcell = 1:nCells(1)
%     for ycell = 1:nCells(2)
        cellCtr = cellCtr+1;
%         ihhi(cellCtr,:,:) = reshape(obj.couplingFilter{ind2sub([nCells(1),nCells(2)],cellCtr)},nCellsTotal,hlen);
    lfilt = length(obj.couplingFilter{xcell}{1});
%     for couplingFilterInd = 1:length(obj.couplingFilter{xcell})
%         ihhi(cellCtr,obj.couplingMatrix{xcell}(couplingFilterInd),:) = obj.couplingFilter{xcell}{couplingFilterInd};
% %             interp1(0:lfilt-1, obj.couplingFilter{xcell}{couplingFilterInd}',[.5+dt:dt:(lfilt)-1], 'linear'); 
%     end
% %     ihhi(cellCtr,cellCtr,:) = obj.postSpikeFilter{xcell};
    ihhi(cellCtr,cellCtr,:) = obj.postSpikeFilter{xcell};%interp1(0:lfilt-1,  obj.postSpikeFilter{xcell},[.5+dt:dt:(lfilt)-1], 'linear'); 
%     end
end
hlen = length(ihhi);
cellCtr = 0;
for xcell = 1:nCells(1)
%     for ycell = 1:nCells(2)
        cellCtr = cellCtr+1;
%         ihhi(cellCtr,:,:) = reshape(obj.couplingFilter{ind2sub([nCells(1),nCells(2)],cellCtr)},nCellsTotal,hlen);
    lfilt = length(obj.couplingFilter{xcell}{1});
    for couplingFilterInd = 1:length(obj.couplingFilter{xcell})
        ihhicpl(cellCtr,obj.couplingMatrix{xcell}(couplingFilterInd),:) = obj.couplingFilter{xcell}{couplingFilterInd};
%             interp1(0:lfilt-1, obj.couplingFilter{xcell}{couplingFilterInd}',[.5+dt:dt:(lfilt)-1], 'linear'); 
    end
% %     ihhi(cellCtr,cellCtr,:) = obj.postSpikeFilter{xcell};
%     ihhi(cellCtr,cellCtr,:) = obj.postSpikeFilter{xcell};%interp1(0:lfilt-1,  obj.postSpikeFilter{xcell},[.5+dt:dt:(lfilt)-1], 'linear'); 
%     end
end

% load('pairspike_1.mat','pairspike');
load('pairspikeall.mat','pairspike');

nCellsTotal = 1;
for xcell = 3%:nCells(1)
for trial = 1:numberTrials
trial
tic
        %%%%%%%%%%%%%%%%%%
        
%         Vstm = log(nlResponse{xcell,ycell});
%         Vstm = (vertcat(obj.nlResponse{:}));
        Vstm = vertcat(obj.linearResponse{:,xcell,1});
        
        
        Vmem = nlfun(interp1([0:slen-1]',Vstm',[.5+dt:dt:slen-1]', 'linear'));
        
        % Vspk = horzcat(spikeResponse{:,1,trial,1});
        
        Ispk = zeros(size(Vmem));
        % Set up simulation dynamics variables
        tsp(1,1:nCellsTotal) = {zeros(round(slen/5),1)};  % allocate space for spike times
        nsp = zeros(1,nCellsTotal);
        jbin = 1;
                
%         tspnext = exprnd(1,1,nCellsTotal);  % Changed for ISETBIO to eliminate toolbox call
        tspnext = rand(1,nCellsTotal);
        rprev = zeros(1,nCellsTotal);
        while jbin <= rlen-hlen;
            tspnext = rand(1,nCellsTotal);
            iinxt = jbin;%jbin:min(jbin+nbinsPerEval-1,rlen);  
            nii = length(iinxt);
            
%             rrnxt = nlfun(Vmem(iinxt,:))*dt/RefreshRate; % Cond Intensity
            for couplingIndex = 1:6
%                 [ispkscpl,jspkscpl] = find(spikeResponse{obj.couplingMatrix{xcell}(couplingIndex),1,trial,1}/dt==jbin);
%                 if ispkscpl ~= 0
%                     Vmem(jbin:jbin+hlen-1,:) = ((Vmem(jbin:jbin+hlen-1,:)).*exp(permute(ihhicpl(xcell,obj.couplingMatrix{xcell}(couplingIndex),1:hlen),[3 2 1])));
%                 end;
                
                ispkscpl = (pairspike{xcell,couplingIndex}(trial,jbin)==1);
                if ispkscpl ~= 0
                    Vmem(jbin:jbin+hlen-1,:) = ((Vmem(jbin:jbin+hlen-1,:)).*exp(permute(ihhicpl(xcell,obj.couplingMatrix{xcell}(couplingIndex),1:hlen),[3 2 1])));
                end;
            end                
            
            rrnxt = (Vmem(iinxt,:))*dt/RefreshRate; % Cond Intensity
            
            rrcum = cumsum(rrnxt+[rprev;zeros(nii-1,nCellsTotal)],1);  % Cumulative intensity
%             if all(tspnext >= rrcum(end,:)) % No spike in this window
%             if all(tspnext <= exp(-rrcum(end,:)))
            rollcomp(jbin,:,trial) = exp(-rrnxt(1,:));
            

            if all(tspnext <= exp(-rrnxt(1,:)))
%             if roll >  exp(-params.bindur*cif(i));  
                jbin = iinxt(end)+1;
                rprev = rrcum(end,:);
            else   % Spike!
%                 [ispks,jspks] =  find(rrcum>=repmat(tspnext,nii,1));
                [ispks,jspks] = find(exp(-rrnxt(1,:))<=(tspnext));
                spcells = unique(jspks(ispks == min(ispks))); % cell number(s)
                ispk = iinxt(min(ispks)); % time bin of spike(s)
                rprev = rrcum(min(ispks),:); % grab accumulated history to here
                
                % Record this spike
                mxi = min(rlen, ispk+hlen); % determine bins for adding h current
                iiPostSpk = ispk+1:mxi;
%                 tspnextall = exprnd(1,1,max(spcells));
                tspnextall = rand(1,max(spcells));
                for ic = 1:length(spcells)
                    
%                     ihhi = reshape(obj.couplingFilter{ind2sub([nCells(1),nCells(2)],ic)},nCellsTotal,hlen);
%                     
                    icell = spcells(ic);
                    nsp(icell) = nsp(icell)+1;
                    tsp{icell}(nsp(icell),1) = ispk*dt;
                    % inhibCell = squeeze(ihhi(:,icell,1:mxi-ispk))';
                    if ~isempty(iiPostSpk)
                        % need to reshape here for when there are no
                        % spikes?
                        
%                         Vmem(iiPostSpk,:) = Vmem(iiPostSpk,:)+permute(ihhi(:,icell,1:mxi-ispk),[3 1 2]);


                        Vmem(iiPostSpk,:) = ((Vmem(iiPostSpk,:)).*exp(permute(ihhi(xcell,xcell,1:mxi-ispk),[3 2 1])));

%                         if nargout == 3  % Record post-spike current
%                             Ispk(iiPostSpk,:)=Ispk(iiPostSpk,:)+permute(ihhi(:,icell,1:mxi-ispk),[3 1 2]);
%                              Ispk(iiPostSpk,:)=Ispk(iiPostSpk,:)+Vmem(iiPostSpk,:);
%                         end
                    end
                    
                    
                    
                    

                
%                 spcells = unique(jspks(ispks == min(ispks))); % cell number(s)
%                 ispk = iinxt(min(ispks)); % time bin of spike(s)
                    
                    rprev(icell) = 0;  % reset this cell's integral
%                     tspnext(icell) = ieExprnd(1,1); % draw RV for next spike in this cell, changed for ISETBIO no toolbox case
                    tspnext(icell) = tspnextall(icell);
%                     tspnext(icell) = rand(1);
%                     VmemAll(ic,:) = Vmem;
                end
                
                %                 [ispks,jspks] =  find(rrcum>=repmat(tspnext,nii,1));

                
                jbin = ispk+1;  % Move to next bin
                % --  Update # of samples per iter ---
                muISI = jbin/(sum(nsp));
                nbinsPerEval = max(120, round(1.5*muISI));
                if isempty(nbinsPerEval); nbinsPerEval = 120; end;
            end
        end
        
        % Remove any extra bins from cell array of spike times
        for j = 1:nCellsTotal
            tsp{j} = tsp{j}(1:nsp(j));
        end
        
        
cellCtr = 0;
for xcell =3%:nCells
%     for ycell = 1:nCells(2)
        
        cellCtr = cellCtr+1;
        spikeTimes{xcell,1,trial,1} = tsp{cellCtr}(1:nsp(cellCtr)); % prune extra zeros
        if size(Vmem,1) > 0
        spikeTimes{xcell,1,trial,2} = Vmem(:,cellCtr); % prune extra zeros
        end
%     end
end
%     end
% end
toc
end%trial
end%xcell
% 


figure; hold on; for ce = 3; for tr = 1:numberTrials; 
if ~isempty(spikeTimes{ce,1,tr,1});
% subplot(6,7,ce); hold on; plot(spikeTimes{ce,1,tr,1},tr,'ok');axis([0 270 0 10]);end;end;end;
subplot(2,1,1);
hold on; line([spikeTimes{ce,1,tr,1},spikeTimes{ce,1,tr,1}],[tr tr-1],'color','k');axis([0 500 0 numberTrials]);end;end;end;

subplot(2,1,2);
convolvewin=gausswin(100);
for trind = 1:numberTrials
    yind= spikeTimes{ce,1,trind,1};
    y(trind,round(yind./dt))=1;
end
% yind=sort(vertcat(spikeTimes{ce,1,1:numberTrials,1}),'ascend');
% y(round(yind./dt))=1;
PSTH_rec=conv(sum(y),convolvewin,'same');
plot(dt:dt:dt*length(PSTH_rec),PSTH_rec);        
axis([0 500 0  max(PSTH_rec)])

% figure; hold on; for ce = 1:1; for tr = 1:10; 
% if ~isempty(spikeTimes{ce,1,tr,1});
% % subplot(6,7,ce); hold on; plot(spikeTimes{ce,1,tr,1},tr,'ok');axis([0 270 0 10]);end;end;end;
% subplot(1,1,ce);
% hold on; line([spikeTimes{ce,1,tr,1},spikeTimes{ce,1,tr,1}],[tr tr-1],'color','k');axis([0 1000 0 10]);end;end;end;
% 
% 
% figure; hold on; for ce = 1:1; for tr = 1:10; 
% if ~isempty(spikeTimes{ce,1,tr,1});
% % subplot(6,7,ce); hold on; plot(spikeTimes{ce,1,tr,1},tr,'ok');axis([0 270 0 10]);end;end;end;
% subplot(1,1,ce);
% hold on; line([spikeResponse{ce,1,tr,1},spikeResponse{ce,1,tr,1}],[tr tr-1],'color','k');axis([0 1000 0 10]);end;end;end;

% 
% figure; hold on; for ce = 1:39; for tr = 1:10; 
% if ~isempty(rgc2.mosaic{1}.spikeResponse{ce,1,tr,1});
% subplot(6,7,ce); hold on; plot(spikeTimes{ce,1,tr,1},tr,'ok');axis([0 270 0 10]);end;end;end;

ph=1;