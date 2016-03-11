function [tsp,Vmem,Ispk] = simGLMcpl(glmprs,Stim);
% [tsp, Vmem,Ispk] = simGLMcpl(glmprs,Stim);
% 
% Compute response of (multi-cell) coupled-glm to stimulus Stim.
%
% Uses time rescaling instead of Bernouli approximation to conditionally
% Poisson process
%
% Dynamics:  Filters the Stimulus with glmprs.k, passes this through a
% nonlinearity to obtain the point-process conditional intensity.  Add a
% post-spike current to the linear input after every spike.

% --------------- Check Inputs ----------------------------------
global RefreshRate;

ncells = size(glmprs.k,3);
nbinsPerEval = 100;  % Default number of bins to update for each spike
dt = glmprs.dt;
% Check that dt evenly divides 1
if mod(1,dt) ~= 0
    dt = 1/round(1/dt);
    fprintf(1, 'glmrunmod: reset dtsim = %.3f\n', dt);
end
if ~isfield(glmprs, 'nlfun')
    nlfun = @exp;
else
    nlfun = glmprs.nlfun;
end

[slen,swid] = size(Stim); % length of stimulus
rlen = round(slen/dt);  % length of binned response

% -------------  Compute filtered resp to signal ----------------
k = glmprs.k;  % stim filter
Vstm = zeros(slen+2,ncells);
if swid == size(k,2) % convolve if k is the same width as stimulus
    for j = 1:ncells
        Vstm(2:end-1,j) = sameconv(Stim,k(:,:,j));
    end
else
    error('Mismatch between stimulus and filter size -- glmrunmod');
end
Vstm = Vstm + repmat(glmprs.dc(:)',slen+2,1);

% -------------- Compute interpolated h current ----------------------
ih = glmprs.ih;
if ~isempty(ih)
    ihthi = [dt:dt:max(glmprs.iht)]';  % time points for sampling
    ihhi = interp1(glmprs.iht, ih, ihthi, 'linear', 0);
    hlen = size(ihhi,1);
    ihhi = permute(ihhi,[1 3 2]); % flip 2nd & 3rd dimensions
else % No ih current
    hlen=1; ihhi=zeros(1,ncells); ihthi = dt;
    if ~isempty(glmprs.iht) 
        warning('No post-spike kernel supplied but iht nonempty!!\n');
    end
end

% -------------  Static nonlinearity & spiking -------------------
Vmem = interp1([0:slen+1]',Vstm,[.5+dt:dt:slen+.5]', 'linear');
if nargout > 2
    Ispk = Vmem*0;
end

% Set up simulation dynamics variables
tsp(1,1:ncells) = {zeros(round(slen/25),1)};  % allocate space for spike times
nsp = zeros(1,ncells);
jbin = 1;

tspnext = exprnd(1,1,ncells);
rprev = zeros(1,ncells);
while jbin <= rlen
    iinxt = jbin:min(jbin+nbinsPerEval-1,rlen);  nii = length(iinxt);
    rrnxt = nlfun(Vmem(iinxt,:))*dt/RefreshRate; % Cond Intensity
    rrcum = cumsum(rrnxt+[rprev;zeros(nii-1,ncells)],1);  % Cumulative intensity
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

% Remove any extra bins from cell array of spike times
for j = 1:ncells
    tsp{j} = tsp{j}(1:nsp(j));
end
