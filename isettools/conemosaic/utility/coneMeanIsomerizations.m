function meanRate = coneMeanIsomerizations(varargin)
% Get the mean photon rate (R*/sec) for the cone types in a mosaic
% 
% Input parameter/value pairs
%   cMosaic
%   pRate    - (x,y,t) of photon absorption rates (R*/sec)
%   coneType - spatial array of LMS cone positions
%
% Return
%   meanRate - 3-vector of R*/sec for each of the cone classes
% 
% 11/2016 JRG (c) isetbio team

%% When we switch to the cMosaic input
%
% p = inputParser; 
% p.addRequired('cMosaic', @(x) isa(x, 'coneMosaic'));
% p.parse(cMosaic,varargin{:});
% cMosaic = p.Results.cMosaic;
% 
% coneType = cMosaic.pattern;
% pRate    = cMosaic.absorptions;

%% Parse parameters
p = inputParser; p.KeepUnmatched = true;
p.addParameter('cMosaic',  [],@(x) isa(x, 'coneMosaic'));
p.addParameter('pRate',    [], @(x) ndims(x) == 3);
p.addParameter('coneType', [], @ismatrix);
p.parse(varargin{:});

cMosaic = p.Results.cMosaic;
pRate = p.Results.pRate;
coneType = p.Results.coneType;

if ~isempty(cMosaic) 
    pRate = cMosaic.absorptions./cMosaic.integrationTime;
    coneType = cMosaic.pattern;
end

%% Get where cones of each type are located in the mosaic
lConeIndices = find(coneType == 2);
mConeIndices = find(coneType == 3);
sConeIndices = find(coneType == 4);

% Reshape (x,y,t) to (x*y,t)
pRateXW = RGB2XWFormat(pRate);
lConeAbsorptions = pRateXW(lConeIndices,:); %#ok<FNDSB>
mConeAbsorptions = pRateXW(mConeIndices,:); %#ok<FNDSB>
sConeAbsorptions = pRateXW(sConeIndices,:); %#ok<FNDSB>

%% Compute means
if ~isempty(lConeAbsorptions), lMean = mean(lConeAbsorptions(:));  
else                           lMean = 0;
end;

if ~isempty(mConeAbsorptions), mMean = mean(mConeAbsorptions(:)); 
else                           mMean = 0;
end;

if ~isempty(sConeAbsorptions), sMean = mean(sConeAbsorptions(:)); 
else                           sMean = 0;
end

meanRate = [lMean, mMean, sMean];

end
