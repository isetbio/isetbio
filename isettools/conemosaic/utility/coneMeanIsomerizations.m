function [lConeMean, mConeMean, sConeMean] = coneMeanIsomerizations(cMosaic, varargin)
% 
% Get the mean photon rate by cone type.
% 
% 11/2016 JRG (c) isetbio team


% parse inputs
p = inputParser; 
p.addRequired('cMosaic', @(x) isa(x, 'coneMosaic'));
p.parse(cMosaic,varargin{:});
cMosaic = p.Results.cMosaic;

coneType = cMosaic.pattern;
pRate    = cMosaic.absorptions;

% Get where cones of each type are located in the mosaic
lConeIndices = find(coneType == 2);
mConeIndices = find(coneType == 3);
sConeIndices = find(coneType == 4);

% Reshape (x,y,t) to (x*y,t)
pRateXW = RGB2XWFormat(pRate);
lConeAbsorptions = pRateXW(lConeIndices,:); %#ok<FNDSB>
mConeAbsorptions = pRateXW(mConeIndices,:); %#ok<FNDSB>
sConeAbsorptions = pRateXW(sConeIndices,:); %#ok<FNDSB>

if ~isempty(lConeAbsorptions); 
    lConeMean = mean(lConeAbsorptions(:)); 
else 
    lConeMean = 0;
end;

if ~isempty(mConeAbsorptions); 
    mConeMean = mean(mConeAbsorptions(:)); 
else
    mConeMean = 0;
end;

if ~isempty(sConeAbsorptions); 
    sConeMean = mean(sConeAbsorptions(:)); 
else
    sConeMean = 0;
end;
