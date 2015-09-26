function spikeTimes = computeSpikesGLM(obj, sensor, outersegment, varargin)
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

nCells = size(obj.spatialRFArray);
nCellsTotal = nCells(1)*nCells(2);
spikeTimes = cell(nCells);

Vstm = obj.nlResponse{1,1};
slen = length(Vstm);
dt = .01; % sensorGet(sensor,'integration time');
% rlen = length([.5+dt:dt:slen+.5]');
rlen = length([.5+dt:dt:slen-1]');
hlen = length(ihcpl{1,1});
nbinsPerEval = 100;
nlfun = @exp;
RefreshRate = 100;

% ihthi = [dt:dt:max(glmprs.iht)]';  % time points for sampling
% ihthi = [dt:dt:slen*dt];
% ihhi = interp1(.001:.01:slen*dt+.01, ih, ihthi, 'linear', 0);
% hlen = length(ihhi);

cellCtr = 0;
for xcell = 1:nCells(1)
    for ycell = 1:nCells(2)
        cellCtr = cellCtr+1;
        ihhi(cellCtr,:,:) = reshape(obj.couplingFilter{ind2sub([nCells(1),nCells(2)],cellCtr)},nCellsTotal,hlen);
  
    end
end

ihhi = reshape(ihhi, hlen, nCellsTotal, nCellsTotal);

        %%%%%%%%%%%%%%%%%%
        
%         Vstm = log(nlResponse{xcell,ycell});
        Vstm = log(vertcat(obj.nlResponse{:}));
        
        
        Vmem = interp1([0:slen-1]',Vstm',[.5+dt:dt:slen-1]', 'linear');
        
        % Set up simulation dynamics variables
        tsp(1,1:nCellsTotal) = {zeros(round(slen/25),1)};  % allocate space for spike times
        nsp = zeros(1,nCellsTotal);
        jbin = 1;
        
        tspnext = exprnd(1,1,nCellsTotal);
        rprev = zeros(1,nCellsTotal);
        
        tspnext = exprnd(1,1,nCellsTotal);
        rprev = zeros(1,nCellsTotal);
        while jbin <= rlen
            iinxt = jbin:min(jbin+nbinsPerEval-1,rlen);  nii = length(iinxt);
            
            rrnxt = nlfun(Vmem(iinxt,:))*dt/RefreshRate; % Cond Intensity
            
%             rrnxt = (Vmem(iinxt))*dt/RefreshRate; % Cond Intensity
            
            rrcum = cumsum(rrnxt+[rprev;zeros(nii-1,nCellsTotal)],1);  % Cumulative intensity
            if all(tspnext >= rrcum(end,:)) % No spike in this window
                jbin = iinxt(end)+1;
                rprev = rrcum(end,:);
            else   % Spike!
                [ispks,jspks] =  find(rrcum>=repmat(tspnext,nii,1));
                spcells = unique(jspks(ispks == min(ispks))); % cell number(s)
                ispk = iinxt(min(ispks)); % time bin of spike(s)
                rprev = rrcum(min(ispks),:); % grab accumulated history to here
                
                % Record this spike
                mxi = min(rlen, ispk+hlen); % determine bins for adding h current
                iiPostSpk = ispk+1:mxi;
                for ic = 1:length(spcells)
                    
%                     ihhi = reshape(obj.couplingFilter{ind2sub([nCells(1),nCells(2)],ic)},nCellsTotal,hlen);
%                     
                    icell = spcells(ic);
                    nsp(icell) = nsp(icell)+1;
                    tsp{icell}(nsp(icell),1) = ispk*dt;
                    if ~isempty(iiPostSpk)
                        Vmem(iiPostSpk,:) = Vmem(iiPostSpk,:)+ihhi(1:mxi-ispk,:,icell);
                        if nargout == 3  % Record post-spike current
                            Ispk(iiPostSpk,:)=Ispk(iiPostSpk,:)+ihhi(1:mxi-ispk,:,icell);
                        end
                    end
                    rprev(icell) = 0;  % reset this cell's integral
                    tspnext(icell) = exprnd(1); % draw RV for next spike in this cell
                end
                jbin = ispk+1;  % Move to next bin
                % --  Update # of samples per iter ---
                muISI = jbin/(sum(nsp));
                nbinsPerEval = max(20, round(1.5*muISI));
            end
        end
        
%         % Remove any extra bins from cell array of spike times
%         for j = 1:nCellsTotal
%             tsp{j} = tsp{j}(1:nsp(j));
%         end
        
        
cellCtr = 0;
for xcell = 1:nCells(1)
    for ycell = 1:nCells(2)
        
        cellCtr = cellCtr+1;
        spikeTimes{xcell,ycell} = tsp{cellCtr}(1:nsp(cellCtr)); % prune extra zeros
        
    end
end
%     end
% end
        ph = 1;