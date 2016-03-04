function spikeTimes = computeSpikesPSF(mosaic, varargin)
% Converts the linear response to % probabilistic spiking output and 
% incorporates a post-spike filter.
%
% Pillow, Shlens, Paninski, Sher, Litke, Chichilnisky, Simoncelli, Nature,
% 2008, licensed for modification, which can be found at
%
% http://pillowlab.princeton.edu/code_GLM.html
% 
% Interoplates the linear response to a finer time scale, converts this 
% interpolated signal to the nonlinear output using the generator function,
% and the nonlinear output to a probabilistic spiking output based on a 
% Poisson process with a refractory period.
% 
% Inputs: the inner retina object with a precomputed linear response
%
% Outputs: the inner retina object with a spiking response
%
% Example:
%   responseSpikes = computeSpikesPSF(ir.mosaic{ii});
% 
% (c) isetbio
% 09/2015 JRG

%%%%%% WRITTEN BY J PILLOW

% -------------  Static nonlinearity & spiking -------------------
% 

nkt = 20;    % Number of time bins in filter;
dt = .01; % Bin size for simulating model & computing likelihood (in units of stimulus frames)

% --- Make basis for post-spike (h) current ------
ihbasprs.ncols = 5;  % Number of basis vectors for post-spike kernel
ihbasprs.hpeaks = [.1 2];  % Peak location for first and last vectors
ihbasprs.b = .5;  % How nonlinear to make spacings
ihbasprs.absref = .1; % absolute refractory period 
[iht,ihbas,ihbasis] = makeBasis_PostSpike(ihbasprs,dt);
psf = ihbasis*[-10 -5 0 2 -2]';  % h current

% Initialize
spResponseSize = size(mosaic.responseLinear{1,1}(:,:,1));
nSamples = size(mosaic.responseLinear{1,1},3);

nCells = size(mosaic.responseLinear);
spikeTimes = cell(nCells);

Vstm = mosaic.responseLinear{1,1};
slen = length(Vstm);

% ihhi = zeros(length(psf),1);
ihhi = psf;

% rlen = length([.5+dt:dt:slen+.5]');
rlen = length([.5+dt:dt:slen-1]');

% hlen = length(ih);
hlen = length(psf);

nbinsPerEval = 100;
nlfun = mosaic.generatorFunction;
RefreshRate = 100;

numberTrials = mosaicGet(mosaic, 'numberTrials');

for trial = 1:numberTrials
    cellCtr = 0;
for xcell = 1:nCells(1)
    for ycell = 1:nCells(2)
        
        Vstm = (mosaic.responseLinear{xcell,ycell});
        % Vstm = vertcat(obj.mosaic.responseLinear{xcell,ycell,1});
        
        nsp = 0;
        tsp = zeros(round(slen/25),1);  % allocate space for spike times
        % Vmem = interp1([0:slen+1]',Vstm,[.5+dt:dt:slen+.5]', 'linear');
        % Vmem = interp1([1:slen],Vstm,[.5+dt:dt:slen+.5], 'linear');
        Vmem = interp1([0:slen-1]',Vstm',[.5+dt:dt:slen-1]', 'linear');

        if (nargout > 2), Ispk = Vmem*0;
        end
        
        jbin = 1; % current time bin
        nsp = 0; % number of spikes
        tspnext = ieExprnd(1,1);  % time of next spike (in rescaled time)
        rprev = 0;  % Integrated rescaled time up to current point
        while jbin <= rlen
            iinxt = jbin:min(jbin+nbinsPerEval-1,rlen);
            % rrnxt = nlfun(Vmem(iinxt))*dt/RefreshRate; % Cond Intensity
            % rrnxt = (Vmem(iinxt))*dt/RefreshRate; % Cond Intensity            
            rrnxt = nlfun(Vmem(iinxt,:))*dt/RefreshRate; % Cond Intensity

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
                tspnext = ieExprnd(1,1);  % draw next spike time
                rprev = 0; % reset integrated intensity
                jbin = ispk+1;  % Move to next bin
                % --  Update # of samples per iter ---
                muISI = jbin/nsp;
                nbinsPerEval = max(20, round(1.5*muISI));
            end
        end
               
        cellCtr = cellCtr+1;
        % Note: x and y indices flipped from normal to match imagesc 
        % This does not happen in computeSpikesGLM because of vertcat on
        % line 58.
        spikeTimes{ycell,xcell,trial,1} = tsp(1:nsp); % prune extra zeros
        if size(Vmem,1) > 0
            spikeTimes{ycell,xcell,trial,2} = Vmem; % prune extra zeros
        end
    end%ycell
end%xcell
end%trial