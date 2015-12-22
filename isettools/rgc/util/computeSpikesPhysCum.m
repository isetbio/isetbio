function spikeTimes = computeSpikesPhys(obj, varargin)
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
hlen = length(ihcpl{1}{1});
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

for trial = 1:numberTrials
trial
tic
cellCtr = 0;

for xcell = 1:nCells(1)
%     for ycell = 1:nCells(2)
        cellCtr = cellCtr+1;
%         ihhi(cellCtr,:,:) = reshape(obj.couplingFilter{ind2sub([nCells(1),nCells(2)],cellCtr)},nCellsTotal,hlen);
    for couplingFilterInd = 1:length(obj.couplingFilter{xcell})
        ihhi(cellCtr,obj.couplingMatrix{xcell}(couplingFilterInd),:) = obj.couplingFilter{xcell}{couplingFilterInd}; 
    end
    ihhi(cellCtr,cellCtr,:) = obj.postSpikeFilter{xcell};
%     end
end

        %%%%%%%%%%%%%%%%%%
        
%         Vstm = log(nlResponse{xcell,ycell});
%         Vstm = (vertcat(obj.nlResponse{:}));
        Vstm = vertcat(obj.linearResponse{:,:,1});
        
        
        Vmem = (interp1([0:slen-1]',Vstm',[.5+dt:dt:slen-1]', 'linear'));
        
        Ispk = zeros(size(Vmem));
        % Set up simulation dynamics variables
        tsp(1,1:nCellsTotal) = {zeros(round(slen/5),1)};  % allocate space for spike times
        nsp = zeros(1,nCellsTotal);
        jbin = 1;
                
        % tspnext = exprnd(1,1,nCellsTotal);  % Changed for ISETBIO to eliminate toolbox call
        tspnext = rand(1,nCellsTotal);
        rprev = zeros(1,nCellsTotal);
        while jbin <= rlen
            iinxt = jbin:min(jbin+nbinsPerEval-1,rlen);  nii = length(iinxt);
            
%             rrnxt = nlfun(Vmem(iinxt,:))*dt/RefreshRate; % Cond Intensity
            
            rrnxt = nlfun(Vmem(iinxt,:))*dt/RefreshRate; % Cond Intensity
            
            rrcum = (cumsum(rrnxt+[rprev;zeros(nii-1,nCellsTotal)],1));  % Cumulative intensity
            % if all(tspnext >= rrcum(end,:)) % No spike in this window
                
            if all(tspnext <= exp(-rrcum(end,:)))
                jbin = iinxt(end)+1;
                rprev = rrcum(end,:);
            else   % Spike!
                rollcomp(jbin,:) = exp(-rrcum(end,:));                 
                [ispks,jspks] =  find(exp(-rrcum)<=repmat(tspnext,nii,1));
                % [ispks,jspks] = find(exp(-rrnxt(1,:))<=(tspnext));
                spcells = unique(jspks(ispks == min(ispks))); % cell number(s)
                ispk = iinxt(min(ispks)); % time bin of spike(s)
                rprev = rrcum(min(ispks),:); % grab accumulated history to here
                
                % Record this spike
                mxi = min(rlen, ispk+hlen); % determine bins for adding h current
                iiPostSpk = ispk+1:mxi;
                % tspnextall = exprnd(1,1,max(spcells));
                tspnextall = rand(1,1,max(spcells));
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
                        % Vmem(iiPostSpk,:) = Vmem(iiPostSpk,:)+permute(ihhi(icell,:,1:mxi-ispk),[3 2 1]);
                        Vmem(iiPostSpk,:) = log(nlfun(Vmem(iiPostSpk,:)).*exp(permute(ihhi(icell,:,1:mxi-ispk),[3 2 1])));
%                         if nargout == 3  % Record post-spike current
%                             Ispk(iiPostSpk,:)=Ispk(iiPostSpk,:)+permute(ihhi(:,icell,1:mxi-ispk),[3 1 2]);
%                              Ispk(iiPostSpk,:)=Ispk(iiPostSpk,:)+Vmem(iiPostSpk,:);
%                         end
                    end
                    rprev(icell) = 0;  % reset this cell's integral
%                     tspnext(icell) = ieExprnd(1,1); % draw RV for next spike in this cell, changed for ISETBIO no toolbox case
                    tspnext(icell) = tspnextall(icell);
%                     VmemAll(ic,:) = Vmem;
                end
                jbin = ispk+1;  % Move to next bin
                % --  Update # of samples per iter ---
                muISI = jbin/(sum(nsp));
                nbinsPerEval = max(20, round(1.5*muISI));
                if isempty(nbinsPerEval); nbinsPerEval = 20; end;
            end
        end
        
        % Remove any extra bins from cell array of spike times
        for j = 1:nCellsTotal
            tsp{j} = tsp{j}(1:nsp(j));
        end
        
        
cellCtr = 0;
for xcell = 1:nCells
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
ph=1;