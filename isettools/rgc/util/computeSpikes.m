function spikeTimes = computeSpikes(nlResponse, ih, sensor, outersegment, varargin)
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



spResponseSize = size(nlResponse{1,1}(:,:,1));
nSamples = size(nlResponse{1,1},3);

nCells = size(nlResponse);
spikeTimes = cell(nCells);

Vstm = nlResponse{1,1};
slen = length(Vstm);
dt = .01; % sensorGet(sensor,'integration time');
% rlen = length([.5+dt:dt:slen+.5]');
rlen = length([.5+dt:dt:slen-1]');
hlen = length(ih);
nbinsPerEval = 100;
nlfun = @exp;
RefreshRate = 100;

% ihthi = [dt:dt:max(glmprs.iht)]';  % time points for sampling
% ihthi = [dt:dt:slen*dt];
% ihhi = interp1(.001:.01:slen*dt+.01, ih, ihthi, 'linear', 0);
% hlen = length(ihhi);
ihhi = ih;
cellCtr = 0;

for trial = 1:10

for xcell = 1:nCells(1)
    for ycell = 1:nCells(2)
        cellCtr = cellCtr + 1;
        Vstm = nlResponse{xcell,ycell};
        
        nsp = 0;
        tsp = zeros(round(slen/25),1);  % allocate space for spike times
        % Vmem = interp1([0:slen+1]',Vstm,[.5+dt:dt:slen+.5]', 'linear');
        % Vmem = interp1([1:slen],Vstm,[.5+dt:dt:slen+.5], 'linear');
        Vmem = interp1([0:slen-1]',Vstm',[.5+dt:dt:slen-1]', 'linear');

        if (nargout > 2), Ispk = Vmem*0;
        end
        
        jbin = 1; % current time bin
        nsp = 0; % number of spikes
        tspnext = exprnd(1);  % time of next spike (in rescaled time)
        rprev = 0;  % Integrated rescaled time up to current point
        while jbin <= rlen
            iinxt = jbin:min(jbin+nbinsPerEval-1,rlen);
            % rrnxt = nlfun(Vmem(iinxt))*dt/RefreshRate; % Cond Intensity
            
            rrnxt = (Vmem(iinxt))*dt/RefreshRate; % Cond Intensity
            rrcum = cumsum(rrnxt)+rprev; % integrated cond intensity
            if (tspnext >= rrcum(end)) % No spike in this window
                jbin = iinxt(end)+1;
                rprev = rrcum(end);
            else   % Spike!
                ispk = iinxt(min(find(rrcum>=tspnext))); % time bin where spike occurred
                nsp = nsp+1;
                tsp(nsp) = ispk*dt; % spike time
                mxi = min(rlen, ispk+hlen); % max time affected by post-spike kernel
                iiPostSpk = ispk+1:mxi; % time bins affected by post-spike kernel
                if ~isempty(iiPostSpk)
                    Vmem(iiPostSpk) = Vmem(iiPostSpk)+ihhi(1:mxi-ispk); %POST SPIKE FIX
                    if nargout == 3  % Record post-spike current
                        Ispk(iiPostSpk) = Ispk(iiPostSpk)+ihhi(1:mxi-ispk);
                    end
                end
                tspnext = exprnd(1);  % draw next spike time
                rprev = 0; % reset integrated intensity
                jbin = ispk+1;  % Move to next bin
                % --  Update # of samples per iter ---
                muISI = jbin/nsp;
                nbinsPerEval = max(20, round(1.5*muISI));
            end
        end
        spikeTimes{xcell,ycell,trial} = tsp(1:nsp); % prune extra zeros
%         if size(Vmem,1) > 0
        spikeTimes{xcell,ycell,trial,2} = Vmem; % prune extra zeros
%         end
        % obj.spkResponse =
        
        ph = 1;
        
    end
end
        ph = 1;
        
end

