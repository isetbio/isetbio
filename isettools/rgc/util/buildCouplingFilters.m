function [couplingFilters, weightMatrix] = buildCouplingFilters(obj, dt)


%%%% Written by J. Pillow

% DTsim = .01; % Bin size for simulating model & computing likelihood (in units of stimulus frames)
nkt = 20;    % Number of time bins in filter;

% dt = .01
% --- Make basis for post-spike (h) current ------
ihbasprs.ncols = 5;  % Number of basis vectors for post-spike kernel
ihbasprs.hpeaks = [.1 2];  % Peak location for first and last vectors
ihbasprs.b = .5;  % How nonlinear to make spacings
ihbasprs.absref = .1; % absolute refractory period 
[iht,ihbas,ihbasis] = makeBasis_PostSpike(ihbasprs,dt);
psf = ihbasis*[-10 -5 0 2 -2]';  % h current

% psf = mosaicGet(obj, 'postSpikeFilter');

% Make some coupling kernels
ihbasprs.ncols = 5;  % Number of basis vectors for post-spike kernel
ihbasprs.hpeaks = [.1 2];  % Peak location for first and last vectors
ihbasprs.b = .5;  % How nonlinear to make spacings
ihbasprs.absref = .1; % absolute refractory period 

[iht,ihbas,ihbasis] = makeBasis_PostSpike(ihbasprs,dt);
hhcpl = ihbasis*[.6;.47;.25;0;0]*2;
hhcpl(:,2) = ihbasis*[-1;-1;0;0;.25]*2;
% ih(:,2,1) = hhcpl(:,2); % 2nd cell coupling to first
% ih(:,1,2) = hhcpl(:,1); % 1st cell coupling to second
ihind = 2;



%%%%

otherCenterLocations = vertcat(obj.cellLocation{:});

nCells = size(obj.cellLocation);

% for mosaic = 1:5
for xcell = 1:nCells(1)
    for ycell = 1:nCells(2)
        
        cellCenterLocation = obj.cellLocation{xcell,ycell};
        
        weightMatrixT = reshape(exp(-(1./(2*obj.rfDiameter))*(sum((bsxfun(@minus, cellCenterLocation, otherCenterLocations)).^2,2).^(1/2))),nCells(1),nCells(2));
        weightMatrixT(weightMatrixT==1) = 0; 
        weightMatrixT(weightMatrixT<0.1) = 0;
        
        
        % weightMatrix{xcell,ycell,mosaic} = weightMatrixT;
        weightMatrix{xcell,ycell} = weightMatrixT;
                
        % couplingFilters{xcell,ycell,mosaic} = reshape(hhcpl(:,1)*weightMatrix(:)',nCells(1),nCells(2),length(hhcpl(:,1)));
        couplingFilters{xcell,ycell} = reshape((hhcpl(:,ihind)*weightMatrixT(:)')',nCells(1),nCells(2),length(hhcpl(:,1)));

        couplingFilters{xcell,ycell}(xcell, ycell, :) = psf;
        
    end
end
% end

% obj = mosaicSet(obj,'couplingFilter', couplingFilters);
% obj = mosaicSet(obj,'couplingMatrix', weightMatrix);